/*
 * Various ugly hacks.  If you're smart, you'll set the compiler
 * options such that this code never gets called.  -JS
 */

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

#include "td3pt.h"

#define DT 10.0e-12

/* Use trapezoid rule for integrals with the DT defined above */

int js_integration_hack(gsl_function *funcp, double T1, double T2, double *Y)
{
  double (*F)(double, void *) = funcp->function;
  void *p                     = funcp->params;
  double x, y;
  int    n, N;

  /* Short-cut */
  if (T2 <= T1)
    return 0.0;

  N = (int)((T2 - T1)/DT + 1e-8);

  if ((T1 + N*DT)/T2 > 1.001 || (T1 + N*DT)/T2 < 0.999)
    {
      mexPrintf("js_integration_hack error 6: T1=%e T2=%e N=%d\n", T1, T2, N);
      return 6;
    }

  y = 0;

  for (n = 0; n <= N; n++)
    {
      x = F(T1 + n*DT, p);
#if (DEBUG)
      mexPrintf("F(%e) = %e\n", n*DT, x);
#endif

      if ((n == 0) || (n == N))
        y +=  x;
      else
        y += 2*x;
    }

  *Y = y * DT/2;

  return 0;
}

