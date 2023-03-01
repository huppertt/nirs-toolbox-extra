/* ----------------------------------------------------------------------
 * Integrand, called by quadpack library.
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

extern double tdgf(const double *, const double *,
                   double, const struct SimParam *);

double tdintegrand(double dt, void *p)
{
  const struct SimParam *sp = (struct SimParam *)p;
  double rs[3], rd[3], phi_src, phi_img;
  int i;

  if (p == NULL)
    mexErrMsgTxt("tdintegrand(): paramter p is NULL");

  /* Make local (modifiable) copies of the src/detector positions */

  for (i = 0; i < 3; i++)
    {
      rs[i] = sp->rsrc[i];
      rd[i] = sp->rdet[i];
    }

  /* The boundary is implicitly at zero.  If slabz is defined, we use
   * its sign to decide which side of the boundary is zero, otherwise
   * we look at the detectors.  If the detectors are also at zero,
   * assume the positive direction is tissue. 
   */

  if      ((sp->slabz < 0) || ((sp->slabz == 0) && (rd[2] < 0)))
    {
      if (rs[2] > 0) rs[2] = 0;    /* Move to interface */
      if (rd[2] > 0) rd[2] = 0;

      rs[2] -= 1/sp->musp;
      rd[2] -= 1/sp->musp;
    }
  else if ((sp->slabz > 0) || ((sp->slabz == 0) && (rd[2] > 0)))
    {
      if (rs[2] < 0) rs[2] = 0;
      if (rd[2] < 0) rd[2] = 0;

      rs[2] += 1/sp->musp;
      rd[2] += 1/sp->musp;
    }
  else
    {
      /* sources and detctors on z=0 plane, assume +Z */

      if (sp->debug)
        mexWarnMsgTxt("Assuming medium is in +Z half-space");

      rs[2] += 1/sp->musp; 
      rd[2] += 1/sp->musp; 
    }

  /* Fluence due to source charge */

  phi_src = tdgf(rs, rd, dt, sp);

  /* Field due to image charge */

  /* (sp->slabz < 0) makes this agree with the matlab code, but I'm
     not convinced it's actually *right* yet */
  if (sp->slabz < 0)
    rs[2] = GET_IMAGE(rs[2], 0 - sp->zext);
  else
    rs[2] = GET_IMAGE(rs[2], 0 + sp->zext);

  phi_img = tdgf(rs, rd, dt, sp);

  /* Turn photon densities into fluences */

  return (sp->v * (phi_src - phi_img));
}
