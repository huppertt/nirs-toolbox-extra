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

#ifndef __TD3PT_H
#define __TD3PT_H

#include "td2pt.h"
#include "gsl/gsl_math.h"

/* Analytic expression for the 3-pt Green's function or calculate by
   integrating up the 2-p Green's functions? */

#define DO_ANALYTIC 1

#define PHI_TOL RELERR/10

enum FLAG { DOMUA, DOMUS };

/* Common short-hand expressions */
/* vector dot product and vector magnitude */

#define DOT(x,y)              (x[0]*y[0] + x[1]*y[1] + x[2]*y[2])
#define NORM(x)               sqrt(DOT(x,x))

/* Subtract, then take the norm of two vectors */

#define SUBNORM(x,y)          sqrt(SQR(x[0] - y[0]) + \
                                   SQR(x[1] - y[1]) + SQR(x[2] - y[2]))

/* Subtract two pairs of vectors, then take their dot product */

#define SUBDOT(x1,x2,y1,y2)  ((x2[0]-x1[0])*(y2[0]-y1[0]) + \
                              (x2[1]-x1[1])*(y2[1]-y1[1]) + \
                              (x2[2]-x1[2])*(y2[2]-y1[2]))

#endif /* __TD3PT_H */
