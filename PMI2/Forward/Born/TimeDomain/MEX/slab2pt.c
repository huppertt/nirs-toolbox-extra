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

#define TOL 1e-5

extern double tdgf(const double *, const double *,
                   double, const struct SimParam *);

double tdintegrand(double dt, void *p)
{
  const struct SimParam *sp = (struct SimParam *)p;
  double phi_src, phi_img, phi_img1, phi_img2;
  double rs[3], rd[3], dphi;
  double z1, zi1, zii1, z2, zi2, zii2;
  int    chg, i;

  if (p == NULL)
    mexErrMsgTxt("tdintegrand(): paramter p is NULL");

  /* Make local (modifiable) copies of the src/detector positions */

  for (i = 0; i < 3; i++)
    {
      rs[i] = sp->rsrc[i];
      rd[i] = sp->rdet[i];
    }

  /* The first physical boundary is implicitly at z=0.  Use the sign
   * of sp->slabz to determine the location of the second physical
   * boundary.
   *
   * Having found the physical boundaries, move the sources one mean
   * free path into the medium.
   */

  if (ABS(sp->slabz) == 0)
    mexErrMsgTxt("Zero slab thickness");
  if (     sp->musp  <= 0)
    mexErrMsgTxt("Non-physical scattering coefficient");

  if (sp->slabz > 0)
    {
      z1 = 0.0;
      z2 = sp->slabz;

      /* Move the sources and detectors to the physical boundary and
         then one mean free path into the medium */

      if (rs[2] < sp->slabz/2)
        rs[2] =         0 + 1/sp->musp;
      else
        rs[2] = sp->slabz - 1/sp->musp;

#if (MOVE_DETS)
      /* Should I be moving detectors here or not?  It's not entirely
         obvious which is the correct approach */

      if (rd[2] < sp->slabz/2)
        rd[2] =         0 + 1/sp->musp;
      else
        rd[2] = sp->slabz - 1/sp->musp;
#endif
    }
  else
    {
      z1 = sp->slabz;
      z2 = 0.0;

      /* Move the sources and detectors to the physical boundary and
         then one mean free path into the medium */

      if (rs[2] > sp->slabz/2)
        rs[2] =         0 - 1/sp->musp;
      else
        rs[2] = sp->slabz + 1/sp->musp;

#if (MOVE_DETS)
      if (rd[2] > sp->slabz/2)
        rd[2] =         0 - 1/sp->musp;
      else
        rd[2] = sp->slabz + 1/sp->musp;
#endif
    }

  /*  Field due to source charge */

  phi_src = tdgf(rs, rd, dt, sp);
  phi_img = 0.0;
  chg     = 1;

  /* For negative times, phi is always zero */

  if (phi_src <= 0.0)
    return 0.0;

  /* Initialize image depths to source position */

  zi1 = rs[2];
  zi2 = rs[2];

  /* Do for a maxiumum of MAXIMG image charges */

  for (i = 1; i < MAXIMG; i++)
    {
      /* Get images of images at opposit boundary */

      /* (sp->slabz < 0) makes this agree with the matlab code, but
         I'm not convinced it's actually *right* yet */
      if (sp->slabz < 0)
        {
          zii1 = GET_IMAGE(zi2, z1 - sp->zext);
          zii2 = GET_IMAGE(zi1, z2 + sp->zext);
        }
      else
        {
          zii1 = GET_IMAGE(zi2, z1 + sp->zext);
          zii2 = GET_IMAGE(zi1, z2 - sp->zext);
        }

      chg  = -chg;

      /* Fluence due to image charges */
      rs[2] = zii1;
      phi_img1 = chg * tdgf(rs, rd, dt, sp);

      rs[2] = zii2;
      phi_img2 = chg * tdgf(rs, rd, dt, sp);

      /* Check for convergence */
      dphi     = ABS(phi_img1 + phi_img2) / ABS(phi_src  + phi_img);
      phi_img += phi_img1 + phi_img2;

      if (dphi < TOL)
        break;

      /* Prepare for next iteration */
      zi1 = zii1;
      zi2 = zii2;
    }

  if (i == MAXIMG)
    mexErrMsgTxt("Fluence failed to converge");

  /* Turn photon densities into fluences */

  return(sp->v * (phi_src + phi_img));
}

