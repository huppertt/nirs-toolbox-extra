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

  if (p == NULL)
    mexErrMsgTxt("tdintegrand(): paramter p is NULL");

  /* tdgf() returns photon density.  I want fluence, multiply by v. */

  return sp->v * tdgf(sp->rsrc, sp->rdet, dt, sp);
}

