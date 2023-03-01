/* ----------------------------------------------------------------------
 *
 * Integrate time-domain two-point Green's function for homogeneous
 * media over the finite width of the detector gate.
 *
 * This code knows nothing about boundary conditions.  Link against
 *  different integrands [td2integrand] to determing which boundary
 *  conditions to use
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
#include "td2pt.h"
#include <gsl/gsl_errno.h>

#define ILIMIT 1024
#define IWORKL (4*ILIMIT)

struct SimParam *sp;

extern void   ierrhnd(const char *, const char *, int, int);
extern void   ierrmsg(int, const char *, double);
extern double tdintegrand(double, void *);
extern double zBnd(void);

#if (!LIBRARY_INTEGRALS)
extern int js_integration_hack(gsl_function *, double, double, double *);
#endif

double TD2GateInt(struct SimParam *spo)
{
  gsl_error_handler_t *olderrh;
  double epsa, epsr, abse=0, phi;
  int    ierr, key;

  /* Make a globally accessible copy of the simulation parameters */

  if ((sp = spo) == NULL)
    mexErrMsgTxt("NULL simulation paramter structure");

  /* Make sure I've replaced the old error handler */

  olderrh = gsl_set_error_handler(ierrhnd);

  /* Integrals with T < 0 are trivial to compute */

  /* Extrapolated boundary length, compute once per call */
  sp->zext = zBnd();

  if (sp->T2 - sp->T1 < 0)
    {
      phi = 0.0;
      ierr = 0;

      if (sp->debug)
	mexPrintf("dt < 0, skipping trivial case");
    }
  else if (sp->T2 == sp->T1)
    {
      /* If the width of the gate is zero, return the value of the
         function instead of its integral */

      phi = tdintegrand(sp->T2, (void *)sp);
      ierr = 0;
    }
  else
    {
      static gsl_integration_workspace *work = NULL;
      gsl_function func;

      if (work == NULL)
	if ((work = gsl_integration_workspace_alloc(IWORKL)) == NULL)
	  mexErrMsgTxt("Error allocating workspace");

      /* Set the integration parameters */
      epsa  = ABSERR;
      epsr  = RELERR;
      abse  = 0.0;

      /* 
       * This is a fairly difficult function to work with, use the
       * most robust algorithm I have available.  
       */

      key = GSL_INTEG_GAUSS51;	  /* 15, 21, 31, 41, 51, or 61 */

      /* 
       * The integrand for t1 < 0 is always going to be zero, that
       * part of the integral is trivial.  
       */

      if (sp->T1 < 0)
	sp->T1 = 0;

      func.function = tdintegrand;
      func.params   = sp;

#if (LIBRARY_INTEGRALS)
      ierr = gsl_integration_qag(&func, sp->T1, sp->T2,
				 epsa, epsr, ILIMIT, key, work, &phi, &abse);
#else
      ierr = js_integration_hack(&func, sp->T1, sp->T2, &phi);
      abse = 1;
#endif
    }

  /* Check for integration errors, print message if needed */

  ierrmsg(ierr, "TD2GateInt", abse);

  /* Restore the original handler */
  gsl_set_error_handler(olderrh);

  if (sp->debug)
    mexPrintf(" integral->%e\n", phi);

  return phi;
}
