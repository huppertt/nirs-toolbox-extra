/* ***********************************************************************
 *
 * Compute two-point Green's functions (Time-Domain DOT)
 *
 * From Matlab:
 *   zBnd = mexzbnd(n, musp, debug);
 *
 * This value is notably more accurate than the one calculated from
 *  inside Matlab.
 *
 * ******************************************************************** */

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

static double calc_zbnd(double n, double musp, int debug);

/* ******************************************************************** */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *N = NULL, *Mu = NULL, *zbnd = NULL;
  int i, nMeas = 0, debug;

  if ((nrhs < 2) || (nrhs > 3))
    mexErrMsgTxt("Incorrect input paramter count.");
  else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments");

  /* Check input types */

  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("First argument [idxRefr] must be a double array");

  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("Second argument [Muspo] must be a double array");

  N  = mxGetPr(prhs[0]);
  Mu = mxGetPr(prhs[1]);

  /* Look for the debug flag */

  debug = 0;

  if (nrhs == 3)
    {
      if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Debug flag must be an integer");
      else
        {
          double *p = mxGetPr(prhs[2]);

          if (p != NULL)
            debug = (p[0] != 0.0);
        }
    }

  if ((mxGetM(prhs[0]) != 1) && (mxGetM(prhs[1]) != 1))
    if (mxGetM(prhs[0]) != mxGetM(prhs[1]))
      mexErrMsgTxt("idxRefr and Musp must be scalar or have same length");
      
  /* Create matrix for the return argument. */

  if (mxGetM(prhs[0]) > mxGetM(prhs[1]))
    nMeas = mxGetM(prhs[0]);
  else
    nMeas = mxGetM(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(nMeas, 1, mxREAL);

  zbnd = mxGetPr(plhs[0]);

  /* Initialize memory */

  for (i = 0; i < nMeas; i++)
    zbnd[i] = 0.0;

  /* Do the actual computation; loop over measurements */

  if (nMeas == 1)
    zbnd[0] = calc_zbnd(N[0], Mu[0], debug);
  else
    {
      for (i = 0; i < nMeas; i++)
	if (mxGetM(prhs[0]) == 1)
	  zbnd[i] = calc_zbnd(N[0], Mu[i], debug);
	else if (mxGetM(prhs[1]) == 1)
	  zbnd[i] = calc_zbnd(N[i], Mu[0], debug);
	else
	  zbnd[i] = calc_zbnd(N[i], Mu[i], debug);
    }

  return;
}

/* ******************************************************************** */

struct SimParam *sp;

static double calc_zbnd(double n, double musp, int debug)
{
  extern double zBnd(void);

  struct SimParam sp0;

  sp = &sp0;

  /* Only these 3 fields are needed */

  sp->n     = n;
  sp->musp  = musp;
  sp->debug = debug;

  return zBnd();
}
