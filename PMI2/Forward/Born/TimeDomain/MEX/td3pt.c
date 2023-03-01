/* ----------------------------------------------------------------------
 *
 * Compute 3-point Green's Function, convolving source and detectors
 *  and integrating over the finite width of the detector gate.
 *
 * This code knows nothing about boundary conditions.  Link against
 *  different integrands [AbsFwd,ScatFwd] to determing which boundary
 *  conditions to use.
 *
 * ------------------------------------------------------------------- */

/* ---------------------------------------------------------------------
 * Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
 *                     Tom Gaudette, Eric Miller, Quan Zhang,
 *                     Jonathan Stott
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * ------------------------------------------------------------------ */

#include <stdlib.h>

#include "td3pt.h"
#include <gsl/gsl_errno.h>

#define ILIMIT 1024
#define IWORKL (4*ILIMIT)

struct SimParam *sp;

extern void   ierrhnd(const char *, const char *, int, int);
extern void   ierrmsg(int, const char *, double);
extern double tdAintegrand(double, void *);
extern double tdSintegrand(double, void *);
extern double  tdAanalytic(double);
extern double  tdSanalytic(double);
extern double zBnd(void);

#if (!LIBRARY_INTEGRALS)
extern int js_integration_hack(gsl_function *, double, double, double *);
#endif

void TDgate(INTEGRAND, double *, const double *, int);
double AbsFwd(double, void *);
double ScatFwd(double, void *);
double SharedFwd(INTEGRAND, double, void *);

/* ------------------------------------------------------------------ */

void TD3GateInt(struct SimParam *spo, 
                  double *A, int nA, const double *rvox, int nvox)
{
  if (nA < nvox)
    mexErrMsgTxt("Forward matrix smaller than nvox");

  if ((nA != nvox) && (nA != 2*nvox))
    mexErrMsgTxt("Forward matrix has inconsistant size");

  /* Make a globally accessible copy of the simulation parameters */

  if ((sp = spo) == NULL)
    mexErrMsgTxt("NULL parameter list");

  sp->zext = zBnd();

  if (sp->muflag[MU_A] && sp->muflag[MU_S])
    {
      TDgate(AbsFwd,  A,      rvox, nvox);
      TDgate(ScatFwd, A+nvox, rvox, nvox);
    }
  else if (sp->muflag[MU_A])
    TDgate(AbsFwd,  A, rvox, nvox);
  else if (sp->muflag[MU_S])
    TDgate(ScatFwd, A, rvox, nvox);

  return;
}

/* ----------------------------------------------------------------------
 * Handle both scattering and absorbing reconstructions, using the
 *  supplied integrand at each voxel
 * ------------------------------------------------------------------- */

void TDgate(INTEGRAND func, double *A, const double *rvox, int nvox)
{
  static gsl_integration_workspace *work = NULL;
  gsl_error_handler_t *olderrh;
  gsl_function ifunc;
  double epsa, epsr, abse;
  int    i, ierr;

  /* Make sure I've replaced the old error handler */
  olderrh = gsl_set_error_handler(ierrhnd);

  /* Only print debugging messages once */

  if (sp->debug && (sp->T1 == sp->T2))
    mexPrintf("T2 = T1, returning integrand\n");
  else if (sp->debug && (sp->T2 - sp->T1 < 0))
    mexPrintf("T2 < T1, trivial integral\n");

#if (LIBRARY_INTEGRALS)
  /* Allocate workspace, if needed */

  if (work == NULL)
    if ((work = gsl_integration_workspace_alloc(ILIMIT)) == NULL)
      mexErrMsgTxt("gsl_integration_workspace_alloc() failed");
#endif

  for (i = 0; i < nvox; i++)
    {
      double phi = 0.0;

      /* Copy current voxel into structure */

      sp->rvox[0] = rvox[3 * i + 0];
      sp->rvox[1] = rvox[3 * i + 1];
      sp->rvox[2] = rvox[3 * i + 2];

      if (sp->T1 == sp->T2)
        {
          /* If the gate width is zero, assume we've been asked for
             the integrand (the instantaneous value at the detector) */

          phi = (*func)(sp->T2, NULL);
        }
      else if (sp->T2 - sp->T1 < 0)
        {
          /* Negative times are trivial to compute as well */

          phi = 0;
        }
      else
        {
          /* For every voxel, integrate over the width of the gate */

          double T1, T2;

          /* Set the integration parameters and integrate */

          ifunc.function = func;
          ifunc.params   = sp;

          epsa = ABSERR;
          epsr = RELERR;

          /* Times < 0 make no contribution */

          T1 = (sp->T1 < 0) ? 0 : sp->T1;
          T2 = (sp->T2 < 0) ? 0 : sp->T2;

#if (LIBRARY_INTEGRALS)
          ierr = gsl_integration_qag(&ifunc, T1, T2, epsa, epsr, 
                                     ILIMIT, GSL_INTEG_GAUSS21, 
                                     work, &phi, &abse);
#else
          ierr = js_integration_hack(&ifunc, T1, T2, &phi);
          abse = 1;
#endif

          /* Check for integration errors, print message if needed */
          ierrmsg(ierr, "TD3GateInt", abse);
        }
      /* Copy results into matrix */
      A[i] = phi;
    }

  /* Restore the original handler */
  gsl_set_error_handler(olderrh);

  return;
}

/************************************************************************
 *
 * Wrapper around the inner integration.  Only the integrands differ.
 *
 ********************************************************************* */

double AbsFwd(double dt, void *p)
{
  /* Increased absorption makes the fluence go down */
  double dmua = -1;
  double phi;

#if (DO_ANALYTIC)
  /* Compute 3-pt Green's function analytically */

  phi = tdAanalytic(dt) * dmua;
#else       /* DO_ANALYTIC */
  /* Compute 3-pt Green's function numerically */

  phi = SharedFwd(tdAintegrand, dt, p) * dmua; 
#endif      /* DO_ANALYTIC */

  return phi;
}

double ScatFwd(double dt, void *p)
{
  const double D = sp->v / (3*(sp->musp + sp->mua));
  double phi;

  /* Actual expression is \grad{\phi}\cdot\grad{G}, not the
     \grad{\phi}\cdot\grad{\phi} computed, so correct for that here.  */

  const double dmus = D / sp->v;

#if (DO_ANALYTIC)
  /* Compute 3-pt Green's function analytically */
  phi = tdSanalytic(dt) * dmus;
#else
  /* Compute 3-pt Green's function numerically */
  phi = SharedFwd(tdSintegrand, dt, p) * dmus;
#endif

  return phi;
}

/* ----------------------------------------------------------------------
 * Common per-voxel code
 * ------------------------------------------------------------------- */

double SharedFwd(INTEGRAND func, double dt, void *p)
{
  static gsl_integration_workspace *work = NULL;
  gsl_function ifunc;
  double epsa, epsr, abse, phi;
  int ierr;

  /* Allocate workspace, if needed */
  if (work == NULL)
    if ((work = gsl_integration_workspace_alloc(ILIMIT)) == NULL)
      mexErrMsgTxt("gsl_integration_workspace_alloc() failed");

  /* Set the base integation paramters and integrate */
      
  epsa = ABSERR;
  epsr = RELERR;
  abse = 0;

  ifunc.function = func;
  ifunc.params   = &dt;

#if (LIBRARY_INTEGRALS)
  ierr = gsl_integration_qags(&ifunc, 0, dt, epsa, epsr,
                              ILIMIT, work, &phi, &abse);
#else
  ierr = js_integration_hack(&ifunc, 0, dt, &phi);
  abse = 1;
#endif

  /* Check for errors */
  ierrmsg(ierr, "SharedFwd()", abse);

  return phi;
}

