/* ----------------------------------------------------------------------
 * Extrapolated boundary for arbitrary index of refraction.  The
 *  extrapolated length is always returned as a positive number.
 *
 * Taken from:
 *    Boundary conditions for the diffusion equation in radiative transfer
 *    Haskell et. al.,  J. Opt. Soc. Am - A   Vol. 11, No. 10,  Oct 1994
 *    pg 2727 - 2741.
 *
 * -------------------------------------------------------------------- */

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

#include <math.h>

#include "td2pt.h"
#include <gsl/gsl_errno.h>

#define ILIMIT 1024                  /* Integration parameters */
#define IWORKL (4*ILIMIT)
#define MAXSUB 256
#define QAGKEY 2

static double nold  = -1.0;
static double rpold =  0.0;
static double rjold =  0.0;

extern struct SimParam *sp;
extern void ierrmsg(int, const char *, double);

static double rfres(double, double);
static double rj_func(double, void *);
static double rphi_func(double, void *);

/* THIS IS DISAGREEING SUBSTANTIALLY (0.155 vs 0.267) WITH
 *  calcExtBnd.m AND I DON'T KNOW WHY!!!!! 
 *
 * CHECK HASKELL's PAPER AND MAKE SURE THE INTEGRALS ARE RIGHT */

double zBnd(void)
{
  double zbnd, Reff, r_phi, r_j, aerr;
  int ierr;

  if (sp == NULL)
    mexErrMsgTxt("zBnd(): extern sp is NULL");

  /*
   * Remember the index of refraction for later.  If it's the same as
   * last invocation, I may be able to lop off LOTS of CPU time by
   * remembering my previous values (this code gets called a lot
   * in the inner loop).
   */

  if ((nold == sp->n) && (rpold > 0) && (rjold > 0))
    {
      r_phi = rpold;
      r_j   = rjold;
    }
  else
    {
      static gsl_integration_workspace *work = NULL;
      gsl_function func;
      double pts[3];
      int npts;

      /* Can't use the old values, recompute from scratch */
      nold = sp->n;

      if (sp->debug)
        mexPrintf("Recomputing zBnd\n");

      /* Compute the critical angle (if one exists) since that's where
       * the discontinuity in the derivative is found.  The integrals
       * on either side of it are fairly easy to integrate. 
       */

      if (sp->n > 1.0)               /* n > 1 has a critical angle */
        {
          pts[0] = 0.0;              /* Lower limit of integration */
          pts[1] = asin(1/sp->n);    /* The critical angle         */
          pts[2] = M_PI/2;           /* Upper limit of integration */
          npts   = 3;
        }
      else
        {
          pts[0] = 0.0;
          pts[1] = M_PI/2;
          npts = 2;
        }

      /* Allocate work space if it hasn't been done already */
      if (work == NULL)
        if ((work = gsl_integration_workspace_alloc(IWORKL)) == NULL)
          mexErrMsgTxt("Error allocating workspace");

      /* Haskell, Eqn 2.3.5(a) */

      func.function = rphi_func;
      func.params   = NULL;

      /* 
       * rphi_func() and rj_func() are a fairly unpleasant function to
       * integrate numerically because while the function is smooth,
       * its derivative is not (the critical angle causes problems).
       * qagp() allows me to subdivide the integral at the
       * discontinuity which gives me a more accurate answer faster
       * than letting the adaptive algorithm find the critical angle.       
       */

      ierr = gsl_integration_qagp(&func, pts, npts, /* 0.0, M_PI/2, */
                                 ABSERR, RELERR, MAXSUB, 
                                 work, &r_phi, &aerr);

      ierrmsg(ierr, "QAGP failed", aerr);

      rpold = r_phi;

      /* Haskell, Eqn 2.3.5(b) */

      func.function = rj_func;
      func.params   = NULL;

      ierr = gsl_integration_qagp(&func, pts, npts, /* 0.0, M_PI/2, */
                                 ABSERR, RELERR, MAXSUB, 
                                 work, &r_j, &aerr);

      ierrmsg(ierr, "QAGP failed", aerr);

      rjold = r_j;
    }

  /* Haskell, eqn 2.3.7 */
  Reff = (r_phi + r_j) / (2 - r_phi + r_j);

  /* Haskell, eqn 2.6.1 */
  zbnd = 2/3.0 * (1 + Reff) / (1 - Reff) / sp->musp;

  /* Always return a positive value; let the caller worry about sign
     conventions */

  return zbnd;
}

/* ----------------------------------------------------------------------
 * Reflection coefficient for photon density
 */

static double rphi_func(double theta, void *dummy)
{
  double costh = cos(theta);
  double sinth = sin(theta);

  return (2 * sinth * costh * rfres(sinth, costh));
}

/* ----------------------------------------------------------------------
 * Reflection coefficient for photon currents
 */

static double rj_func(double theta, void *dummy)
{
  double costh = cos(theta);
  double sinth = sin(theta);

  return (3 * sinth * costh * costh * rfres(sinth, costh));
}

/* ----------------------------------------------------------------------
 * Fresnel coefficients
 */

static double rfres(double sinth, double costh)
{
  double rf1, rf2, r_fres, csthp;

  csthp = 1.0 - SQR(nold * sinth);

  if (csthp <= 0.0)
    r_fres = 1.0;                    /* Total internal reflection */
  else
    {
      csthp = sqrt(csthp);

      rf1   = (nold * csthp - costh) / (nold * csthp + costh);
      rf2   = (nold * costh - csthp) / (nold * costh + csthp);

      r_fres= (SQR(rf1) + SQR(rf2)) / 2;
    }

  return r_fres;
}
