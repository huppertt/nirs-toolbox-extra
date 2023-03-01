/* *********************************************************************
 *
 * Common declarations used by all the different C routines.
 *
 ******************************************************************** */

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

#ifndef _TD2PT_H
#define _TD2PT_H

#include <math.h>
#include "gsl/gsl_integration.h"

#include "mex.h"

/* Use proper integration routines instead of my hacked routine */
#define LIBRARY_INTEGRALS 1

/* Do I move detectors 1 scattering length into the medium or not? */
#define MOVE_DETS 1

/* Vacuum sped of light in cm/s */
#define C_VACUUM 2.99792458e10

/* Maximum allowed number of image charges */
#define MAXIMG 50

/* ******************************************************************** */

#define OFFSET_SRC   (1-1)
#define OFFSET_DET   (2-1)
#define OFFSET_WVL   (4-1)
#define OFFSET_DELAY (6-1)
#define OFFSET_GATE  (7-1)

/* Common short-hand expressions */

#ifndef SQR
#define SQR(x) gsl_pow_2(x)
#endif

#ifndef CUBE
#define CUBE(x) gsl_pow_3(x)
#endif

/* ******************************************************************** */

#define MU_A 0
#define MU_S 1

#define ABSERR 1.0e-4
#define RELERR ABSERR

typedef double (* INTEGRAND)(double, void *);

struct SimParam 
{
  double rsrc[3], rdet[3], rvox[3];
  double T1, T2;
  double musp, mua, g, n, v;
  double slabz, zext;
  int muflag[2], debug;
};

/* The library supplies M_PI if it doesn't already exist */

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#define GET_IMAGE(ZS,ZB) (ZS - 2*(ZS - ZB))

#endif /* _TD2PT_H */

