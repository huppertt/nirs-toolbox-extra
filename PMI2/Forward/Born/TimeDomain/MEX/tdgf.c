/* ----------------------------------------------------------------------
 *
 * Time domain Green's function for DPDW (see Patterson, 1989).
 * NOTE: This definition of D differs from his by a factor of 'c'
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

#include <math.h>
#include "td2pt.h"
#include "gsl/gsl_machine.h"

#define POW32(x)  sqrt((x)*(x)*(x))

double tdgf(const double rsrc[3], const double rdet[3],
	    double dt, const struct SimParam *sp)
{
  /* Time-independant diffusion coefficient */
  const double D = sp->v / (3 * (sp->musp + sp->mua));

  double rsqr, arg;

  if (sp == NULL)
    mexErrMsgTxt("tdgf(): sp is NULL");

  if (dt < 0.0)
    {
      /* Enforce causality */

#if (DEBUG)
      mexPrintf("Non-causal times, returns 0\n");
#endif
      return 0.0;
    }
  else if (dt == 0.0)
    {
      /* Causal, but numerically unstable.  Green's function does go
         to zero though, in the limit t->0, at least analytically for
         r^2 > 0. */
#if (DEBUG)
      mexPrintf("Substituting analytical result for t=0\n");
#endif
      return 0.0;
    }
  else
    {
      /* Calculate S-D separation */
      rsqr = SQR(rdet[0] - rsrc[0]) +
             SQR(rdet[1] - rsrc[1]) + SQR(rdet[2] - rsrc[2]);

      /* Argument of the exponential */
      arg = rsqr / (4 * D * dt) + sp->mua * sp->v * dt;

      if (arg <= 0.0)
        mexErrMsgTxt("Illegal optical properties in tdgf()");

      /* GSL_LOG_DBL_MAX = log(1E+308) or so, anything more than this
       * will be truncated to zero when we do exp(-arg), so save
       * ourselves the call to the exponential and return zero
       * directly*/

      if (arg > GSL_LOG_DBL_MAX)
        {
#if (DEBUG)
          mexPrintf("Range error, rounding to zero\n");
#endif
          return 0.0;
        }
      else
        {
          /* Compute the Green's function (photon density) */
          return (exp(-arg) / POW32(4 * M_PI * D * dt));
        }
    }
}


