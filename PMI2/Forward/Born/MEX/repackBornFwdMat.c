/* ***********************************************************************
 *
 * Go from packed dense representation to correct sparse unpacked
 * representation.  For some reason, Matlab never seems to get the
 * ordering right which means that a fairly simple matrix copy ends up
 * taking forever and a day.
 *
 * From Matlab:
 *   A = repackBornFwdMat(SD, MeasList, A, OptProp, Debug);
 *
 * Jonathan Stott, 2003
 *
 * ******************************************************************** */

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

/* ******************************************************************** */

static mxArray *repack1Matrix(const mxArray *, const mxArray *, 
                              const mxArray *, int debug);
static mxArray *repack2Matrix(const mxArray *, const mxArray *, 
                              const mxArray *, int debug);

/* ******************************************************************** */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *Lambda = NULL, *MeasList = NULL, *A0 = NULL;
  double        *MuVec = NULL;
  int debug = 0;

  if ((nrhs != 4) && (nrhs != 5))
    mexErrMsgTxt("Incorrect input paramter count.");
  else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments");

  /* Check input types */

  if (!mxIsStruct(prhs[0]))
    mexErrMsgTxt("First argument [SD] must be a structure");
  else
    {
      const mxArray *SD = prhs[0];

      if (SD == NULL)
        mexErrMsgTxt("Assertion failed, SD==NULL");

      if ((Lambda = mxGetField(SD, 0, "Lambda")) == NULL)
        mexErrMsgTxt("SD.Lambda not found");
    }

  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("Second argument [MeasList] must be an array");
  else
    {
      MeasList = prhs[1];

      if (MeasList == NULL)
        mexErrMsgTxt("Assertion failed, MeasList==NULL");
    }

  if (mxIsSparse(prhs[2]))
    mexErrMsgTxt("Forward matrix [A] is already sparse");
  else if (!mxIsDouble(prhs[2]))
    mexErrMsgTxt("Third argument [A] must be an array");
  else
    {
      A0 = prhs[2];

      if (A0 == NULL)
        mexErrMsgTxt("Assertion failed, A==NULL");
    }

  if (!mxIsDouble(prhs[3]))
    mexErrMsgTxt("Fourth argument [OptProp] must be an array");
  else
    {
      if ((MuVec = mxGetPr(prhs[3])) == NULL)
        mexErrMsgTxt("Assertion failed, OptProp==NULL");
    }

  /* Look for the debug flag */

  if (nrhs == 5)
    {
      if (!mxIsDouble(prhs[4]))
        mexErrMsgTxt("Debug flag must be an integer");
      else
        {
          double *p = mxGetPr(prhs[4]);

          if (p != NULL)
            debug = (p[0] != 0.0);
        }
    }

  /* Do the actual computation; loops over measurements and voxels*/

  if (MuVec[0] != 0 && MuVec[1] != 0)
    plhs[0] = repack2Matrix(Lambda, MeasList, A0, debug);
  else
    plhs[0] = repack1Matrix(Lambda, MeasList, A0, debug);

  return;
}

/* *********************************************************************
 *
 * Repack a matrix when either calc_mua or calc_mus is true, but not
 * both calc_mua and calc_mus true.
 *
 * ****************************************************************** */

static mxArray *repack1Matrix(const mxArray *Lambda, const mxArray *MeasList, 
                              const mxArray *A0, int debug)
{
  mxArray *Ap;
  double  *Apr, *MLp, *pA0;
  int     nMeas, mlwth, nVoxl, nWvl, *Ajc, *Air;
  int     iw, iv, im, isparse;

  nWvl  = mxGetM(Lambda) * mxGetN(Lambda);

  nMeas = mxGetM(MeasList);
  mlwth = mxGetN(MeasList);

  if (mlwth < 4)
    mexErrMsgTxt("MeasList must be an N x 4+ array");

  if (mxGetM(A0) != nMeas)
    mexErrMsgTxt("Size of MeasList and A do not match");

  nVoxl = mxGetN(A0);

  MLp = mxGetPr(MeasList);
  pA0 = mxGetPr(A0);

  if ((MLp == NULL) || (pA0 == NULL))
    mexErrMsgTxt("Assertion failed, MLp or pA0 NULL");

  if (debug)
    mexPrintf("Block packing forward matrix\n");

  if ((Ap = mxCreateSparse(nMeas, nWvl*nVoxl, 
                           nMeas*nVoxl, mxREAL)) == NULL)
    mexErrMsgTxt("Unable to allocate sparse array");

  if (debug)
    mexPrintf("Allocated sparse storage\n");

  Apr = mxGetPr(Ap);
  Air = mxGetIr(Ap);
  Ajc = mxGetJc(Ap);

  if ((Apr == NULL) || (Air == NULL) || (Ajc == NULL))
    mexErrMsgTxt("Assertion failed, unable to access sparse vectors");

  isparse = 0;

  for (iw = 0; iw < nWvl; iw++)
    {
      for (iv = 0; iv < nVoxl; iv++)
        {
          int spN = iw * nVoxl + iv;
          int AN = iv;

          Ajc[spN] = isparse;

          for (im = 0; im < nMeas; im++)
            {
              /* iWvl = MeasList(im,4); */
              int iWvl = MLp[nMeas * 3 + im] - 1;

              if ((iWvl < 0) || (iWvl >= nWvl))
                mexErrMsgTxt("Wavelength > length(SD.Lambda)");

              if (isparse > nMeas * nVoxl)
                {
                  mexPrintf("nw=%d, nv=%d, nm=%d\n", nWvl, nVoxl, nMeas);
                  mexPrintf("iw=%d, iv=%d, im=%d, isparse=%d\n",
                            iw, iv, im, isparse);
                  mexErrMsgTxt("Fell off end of sparse array");
                }

              if (iWvl == iw)
                {
                  int spM = im;
                  int AM  = im;

                  /* Update the tables */

                  Air[isparse] = spM;
                  Apr[isparse] = pA0[AN * nMeas + AM];

                  isparse += 1;
                }
            }
        }

      if (debug)
        mexPrintf("Wavelength %d\n", iw+1);
    }

  if (isparse != nMeas * nVoxl)
    {
      mexPrintf("isparse = %d vs %d\n", isparse, nMeas * nVoxl);
      mexErrMsgTxt("Element count off");
    }

  /* Tack on final count of elements */

  Ajc[nWvl*nVoxl] = isparse;

  return Ap;
}

/* *********************************************************************
 *
 * Repack a matrix when both calc_mua and calc_mus are true.
 *
 * ****************************************************************** */

static mxArray *repack2Matrix(const mxArray *Lambda, const mxArray *MeasList, 
                              const mxArray *A0, int debug)
{
  mxArray *Ap;
  double  *Apr, *MLp, *pA0;
  int     nMeas, mlwth, nVoxl, nWvl, *Ajc, *Air;
  int     iw, iv, im, isparse;

  nWvl  = mxGetM(Lambda) * mxGetN(Lambda);

  nMeas = mxGetM(MeasList);
  mlwth = mxGetN(MeasList);

  if (mlwth < 4)
    mexErrMsgTxt("MeasList must be an N x 4+ array");

  if (mxGetM(A0) != nMeas)
    mexErrMsgTxt("Size of MeasList and A do not match");

  nVoxl = mxGetN(A0) / 2;

  MLp = mxGetPr(MeasList);
  pA0 = mxGetPr(A0);

  if ((MLp == NULL) || (pA0 == NULL))
    mexErrMsgTxt("Assertion failed, MLp or pA0 NULL");

  if (debug)
    mexPrintf("Block packing forward matrix\n");

  if ((Ap = mxCreateSparse(nMeas, nWvl*2*nVoxl, 
                           nMeas*2*nVoxl, mxREAL)) == NULL)
    mexErrMsgTxt("Unable to allocate sparse array");

  if (debug)
    mexPrintf("Allocated sparse storage\n");

  Apr = mxGetPr(Ap);
  Air = mxGetIr(Ap);
  Ajc = mxGetJc(Ap);

  if ((Apr == NULL) || (Air == NULL) || (Ajc == NULL))
    mexErrMsgTxt("Assertion failed, unable to access sparse vectors");

  isparse = 0;

  for (iw = 0; iw < nWvl; iw++)
    {
      for (iv = 0; iv < 2*nVoxl; iv++)
        {
          int spN = iw * (2*nVoxl) + iv;
          int AN  = iv;

          Ajc[spN] = isparse;

          for (im = 0; im < nMeas; im++)
            {
              /* iWvl = MeasList(im,4); */
              int iWvl = MLp[nMeas * 3 + im] - 1;

              if ((iWvl < 0) || (iWvl >= nWvl))
                mexErrMsgTxt("Wavelength > length(SD.Lambda)");

              if (isparse > nMeas * (2*nVoxl))
                {
                  mexPrintf("nw=%d, nv=%d, nm=%d\n", nWvl, 2*nVoxl, nMeas);
                  mexPrintf("iw=%d, iv=%d, im=%d, isparse=%d\n",
                            iw, iv, im, isparse);
                  mexErrMsgTxt("Fell off end of sparse array");
                }

              if (iWvl == iw)
                {
                  int spM = im;
                  int AM  = im;

                  /* Update the tables */

                  Air[isparse] = spM;
                  Apr[isparse] = pA0[AN * nMeas + AM];

                  isparse += 1;
                }
            }
        }

      if (debug)
        mexPrintf("Wavelength %d\n", iw+1);
    }

  if (isparse != nMeas * 2*nVoxl)
    {
      mexPrintf("isparse = %d vs %d\n", isparse, nMeas * 2*nVoxl);
      mexErrMsgTxt("Element count off");
    }

  /* Tack on final count of elements */

  Ajc[nWvl * 2*nVoxl] = isparse;

  return Ap;
}
