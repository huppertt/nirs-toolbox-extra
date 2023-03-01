/* ---------------------------------------------------------------------
 *
 * Compute 3-point Green's Function, G(rv-rs,t'-0) * G(rd-rv, t-t'), for
 *  scattering and absorbing perturbations.  Infinite medium.
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
extern double tdAAnalytBase(double dt, double rs[3], double rd[3]);
extern double tdSAnalytBase(double dt, double rs[3], double rd[3]);

extern struct SimParam *sp;

/* ----------------------------------------------------------------------
 * Absorbtion Integrand.  Start time implicitly 0, end time passed
 *  as paramter p, intermediate time tvox explicit argument.
 * ------------------------------------------------------------------- */

double tdAintegrand(double tvox, void *p)
{
  double tdet = *(double *)p;
  double phi1, phi2;

  /* Non-positive times are trivial */

  if ((tvox <= 0) || (tdet - tvox <= 0))
    return 0.0;

  /* Call Green's functions.  Green's function returns photon
   densities, I need fluences */

  phi1 = sp->v * tdgf(sp->rsrc, sp->rvox, tvox,        sp);
  phi2 = sp->v * tdgf(sp->rdet, sp->rvox, tdet - tvox, sp);

  return (phi1 * phi2);
}

/* ----------------------------------------------------------------------
 * Scattering Integrand.  Start time implicitly 0, end time passed
 *  as paramter p, intermediate time tvox explicit argument.
 * ------------------------------------------------------------------- */

double tdSintegrand(double tvox, void *p)
{
  double tdet = *(double *)p;
  double phi1[3], phi2[3], *G, grad;

  /* Non-positive times are trivial */

  if ((tvox <= 0) || (tdet - tvox <= 0))
    return 0.0;

  /* Call Green's functions */

  G = tdggf(sp->rsrc, sp->rvox, tvox,        sp);
  memcpy(phi1, G, 3*sizeof(double));

  G = tdggf(sp->rdet, sp->rvox, tdet - tvox, sp);
  memcpy(phi2, G, 3*sizeof(double));

  grad = DOT(phi1, phi2);

  /* Green's function returns photon densities, I need fluences */

  grad = sp->v * grad * sp->v;
  return grad;
}

/* ****************************************************************** */

double tdAanalytic(double dt)
{
  return tdAAnalytBase(dt, sp->rsrc, sp->rdet);
}

double tdSanalytic(double dt)
{
  return tdSAnalytBase(dt, sp->rsrc, sp->rdet);
}
