/* ---------------------------------------------------------------------
 *
 * Compute analytic solution to 3pt Green's function.  Image charges
 *  (if any) and source motion handled by caller.
 *  
 * ------------------------------------------------------------------ */

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
#include <string.h>

#include "td3pt.h"

extern double   tdgf(const double *, const double *,
                     double, const struct SimParam *);
extern double *tdggf(const double *, const double *,
                     double, const struct SimParam *);

extern struct SimParam *sp;

double tdAAnalytBase(double dt, double rs[3], double rd[3])
{
  const double D = sp->v / (3 * (sp->musp + sp->mua));
  double tmps[3], tmpd[3], r1, r2, phi;

  r1 = SUBNORM(sp->rvox, rd);
  r2 = SUBNORM(sp->rvox, rs);

  /* divide by zeros are bad */

  if (r1 <= 0.0 || r2 <= 0.0)
    mexErrMsgTxt("rsrc-rvox or rvox-rdet is zero");

  if (dt <= 0.0)
    {
#if (DEBUG)
      mexWarnMsgTxt("dt <= 0");
#endif
      return 0;
    }

  /* Initialize fake src/det vectors */

  tmps[0] = tmps[1] = tmps[2] = 0.0;
  memcpy(tmpd, tmps, 3*sizeof(double));  tmpd[2] = r1 + r2;

  /* Put it all together */

  phi  = sp->v * tdgf(tmps, tmpd, dt, sp);
  phi *= sp->v * (r1 + r2) / (4*M_PI*D * r1 * r2);

  return phi;
}

double tdSAnalytBase(double dt, double rs[3], double rd[3])
{
  const double D = sp->v / (3 * (sp->musp + sp->mua));
  double r1, r2, phi, gdotg;

  r1 = SUBNORM(sp->rvox, rd);
  r2 = SUBNORM(sp->rvox, rs);

  /* divide by zeros are bad */

  if (r1 <= 0.0 || r2 <= 0.0)
    mexErrMsgTxt("rsrc-rvox or rvox-rdet is zero");

  if (dt <= 0.0)
    {
#if (DEBUG)
      mexWarnMsgTxt("dt <= 0");
#endif
      return 0;
    }

  /* Angular term, dot product of two unit vectors  */

  gdotg = SUBDOT(sp->rvox, rs, sp->rvox, rd) / (r1 * r2);

  /* Put it all together */

  phi  = gdotg * tdAAnalytBase(dt, rs, rd);
  phi *= (SQR(r1) - r1*r2 + SQR(r2)) / (2 * D * dt * r1 * r2)
       + SQR(r1 + r2) / SQR(2 * D * dt);

  return phi;
}
