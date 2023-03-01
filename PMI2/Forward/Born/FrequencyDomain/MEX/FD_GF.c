/* *********************************************************************
 * Try and compute Green's functions a bit faster
 *
 * Phi = mexFDGF(rSrc, rVox, D, K, V);
 *
 * Over-rides the Matlab FD_GF.m routine if placed first in the path,
 *  doing so provides a 4x speed boost with no code changes.
 * ******************************************************************* */

/* ---------------------------------------------------------------------
 * Copyright (C) 2004, David Boas, Dana Brooks, Rick Gaudette, 
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

#include "mex.h"
#include "math.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) <=(y) ? (x) : (y))
#define SQR(x)   ((x)*(x))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *mPhi;
  double D, Kr, Ki, V, *phiR, *phiI, *p;
  double *rsrc, *rvox;
  int    i, ns1=0, ns2, nv1=0, nv2;
  mxComplexity isCplx;

  /* ********************* Parse the inputs ********************* */

  if (nrhs != 5)
    mexErrMsgTxt("Incorrect input parameter count.");
  else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments");

  /* First argument, source position vector */

  if ((rsrc = mxGetPr(prhs[0])) == NULL)
    mexErrMsgTxt("Error extracting source positions");
  else
    {
      ns1 = mxGetM(prhs[0]);
      ns2 = mxGetN(prhs[0]);

      if (ns2 != 3)
        mexErrMsgTxt("rSrc is not an Nx3 positions vector");
    }

  /* Second argument, detector position vector */

  if ((rvox = mxGetPr(prhs[1])) == NULL)
    mexErrMsgTxt("Error extracting voxel positions");
  else
    {
      nv1 = mxGetM(prhs[1]);
      nv2 = mxGetN(prhs[1]);

      if (nv2 != 3)
        mexErrMsgTxt("rVox is not an Nx3 positions vector");
    }

  /* Third argument, diffusion coefficient */

  if ((p = mxGetPr(prhs[2])) == NULL)
    mexErrMsgTxt("Error extracting diffusion coefficient");

  D = p[0];

  if (D <= 0)
    mexErrMsgTxt("Diffusion coefficient must be strictly positive");

  /* Fourth argument, diffuse wave-vector */

  if ((p = mxGetPr(prhs[3])) == NULL)
    mexErrMsgTxt("Error extracting real part of K");

  Kr = p[0];

  if ((p = mxGetPi(prhs[3])) == NULL)
    Ki = 0.0;
  else
    Ki = p[0];

  /* Fifth argument, speed of light */

  if ((p = mxGetPr(prhs[4])) == NULL)
    mexErrMsgTxt("Error extracting velocity");

  V = p[0];

  if (V <= 0)
    mexErrMsgTxt("Velocity must be strictly positive");

  /* If there's no return argument, there's no reason to do the rest
   * of the calculation.  Keep the input checks, though, as that could
   * conceivably be useful by itself. 
   */

  if (nlhs == 0) 
    return;

  /* ****************** Calculate the fluence ******************* */

  isCplx = (Ki == 0) ? mxREAL : mxCOMPLEX;

  mPhi = mxCreateDoubleMatrix(MAX(ns1,nv1), 1, isCplx);
  phiR = mxGetPr(mPhi);
  phiI = mxGetPi(mPhi);

  for (i = 0; i < MAX(ns1, nv1); i++)
    {
      int    iS, iV;
      double dR;

      /* Index into source/detector positions table */

      iS = (ns1 == 1) ? 0 : i;
      iV = (nv1 == 1) ? 0 : i;

      /* Interval at this measurement */

      dR = sqrt(SQR(rsrc[0*ns1 + iS] - rvox[0*nv1 + iV]) +
                SQR(rsrc[1*ns1 + iS] - rvox[1*nv1 + iV]) +
                SQR(rsrc[2*ns1 + iS] - rvox[2*nv1 + iV]));

      /* Calculate the Green's function */

      phiR[i] = exp(-Kr * dR) / (4 * M_PI * D * dR);

      if (phiI != NULL)
        {
          phiI[i] = phiR[i] * sin(-Ki * dR);
          phiR[i] = phiR[i] * cos(-Ki * dR);
        }
    }

  /* ************** Pass back the calcuated values ************** */

  if (nlhs > 0)
    plhs[0] = mPhi;

  return;
}
