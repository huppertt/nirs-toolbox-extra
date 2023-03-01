/* ----------------------------------------------------------------------
 *
 * Gradient of Time domain Green's function for DPDW (see Patterson, 1989).
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

#ifndef POW32
#define POW32(x)  sqrt((x)*(x)*(x))
#endif

extern double tdgf(const double[3], const double[3],
                   double, const struct SimParam *);

double *tdggf(const double rsrc[3], const double rdet[3], 
              double dt, const struct SimParam *sp)
{
  static double G[3];
  double D, rsqr, phi, grad;

  if (sp == NULL)
    mexErrMsgTxt("tdggf(): sp is NULL");

  G[0] = G[1] = G[2] = 0.0;

  if (dt <= 0.0)
    {
      /* Enforce causality */

#if (DEBUG)
      mexPrintf("Non-causal times\n");
#endif
      return G;
    }

  /* S-D separation */
  rsqr = SQR(rdet[0] - rsrc[0]) + 
         SQR(rdet[1] - rsrc[1]) + SQR(rdet[2] - rsrc[2]);

  /* Time-independant diffusion coefficient */
  D = sp->v / (3 * (sp->musp + sp->mua));

  /* Compute the Green's function */
  phi = tdgf(rsrc, rdet, dt, sp);

  /* Compute analytical gradients from basic Green's function */

  grad = phi / (-2 * D * dt);

  G[0] = (rdet[0] - rsrc[0]) * grad;
  G[1] = (rdet[1] - rsrc[1]) * grad;
  G[2] = (rdet[2] - rsrc[2]) * grad;

  return G;
}

