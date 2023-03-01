/* ***********************************************************************
 *
 * Compute three-point Green's functions (Time-Domain DOT)
 *
 * From Matlab:
 *   A = mextd3ptXX(SD, Medium, MeasList, OptProp, Debug);
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
#include <string.h>
#include "td3pt.h"

#ifndef SQR
#define SQR(X) ((X)*(X))
#endif

/* ******************************************************************** */

static void TD3pt(const mxArray *SD, const mxArray *Medium, 
                  const mxArray *ML, const mxArray *OP,
                  const double *rv, int nv, int i, double *A, int debug);

extern void TD3GateInt(struct SimParam *, double *, int, const double *, int);
extern double *getVector(const mxArray *, const char *, int *);
extern double  getField(const mxArray *,  const char *, unsigned int);
extern int     getSampleVolume(const mxArray *, double **, double **);

/* ******************************************************************** */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *SD = NULL, *Medium = NULL, *MeasList = NULL, *OP = NULL;
  int i, iMeas, nMeas = 0, nVox = 0, nAVox, debug = 0;
  double *A, *q, *rvox = NULL, *vvox = NULL, *Atmp;
  mxArray *p, *Atmpm;

  if ((nrhs < 2) || (nrhs > 6))
    mexErrMsgTxt("Incorrect input paramter count.");
  else if (nlhs > 1)
    mexErrMsgTxt("Too many output arguments");

  /* Check input types */

  if (!mxIsStruct(prhs[0]))
    mexErrMsgTxt("First argument [SD] must be a structure");
  else
    SD = prhs[0];

  if (!mxIsStruct(prhs[1]))
    mexErrMsgTxt("Second argument [Medium] must be a structure");
  else
    Medium = prhs[1];

  /* MeasList = NULL; */

  if (nrhs > 2)
    {
      int mlwth;

      if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("Third argument [MeasList] must be an array");

      nMeas = mxGetM(prhs[2]);
      mlwth = mxGetN(prhs[2]);

      if (nMeas > 0)
        {
          if (mlwth < 8)
            mexErrMsgTxt("MeasList must be an N x 8+ array or empty");
          else
            MeasList = prhs[2];
        }
      /* Else leave MeasList NULL and fall through below */
    }

  /* Try the backup copy in SD if I don't have a measurement list yet */

  if (MeasList == NULL)
    {
      int mlwth;

      if ((MeasList = mxGetField(SD, 0, "MeasList")) == NULL)
        mexErrMsgTxt("MeasList not passed and SD.MeasList not found");

      if (!mxIsDouble(MeasList))
        mexErrMsgTxt("SD.MeasList must be an array");

      nMeas = mxGetM(MeasList);
      mlwth = mxGetN(MeasList);

      if (mlwth < 8)
        mexErrMsgTxt("SD.MeasList must be an N x 8+ array");
    }

  if (!mxIsDouble(prhs[3]))
    mexErrMsgTxt("Fourth argument [OptProp] must be an array");
  else
    {
      int nOP;

      OP  = prhs[3];
      nOP = mxGetNumberOfElements(OP);

      if (nOP != 2)
        mexErrMsgTxt("OptProp must be a 2-element vector");
    }

  if (nrhs == 5)
    {
      if (!mxIsDouble(prhs[4]))
        mexErrMsgTxt("Debug flag must be an integer");
      else
        {
          double *p = mxGetPr(prhs[4]);

          debug = (p[0] != 0.0);
        }
    }

  /* Determine size of return matrix and create a vector of voxel positions */

  if ((p = mxGetField(Medium, 0, "CompVol")) == NULL)
    mexErrMsgTxt("Medium.CompVol - field not found");
  else
    nVox = getSampleVolume(p, &rvox, &vvox);
  
  /* Create matrix for the return argument.  The matrix is generated
  *   as the transpose of the actual forward matrix because that's
  *   MUCH easier to work with inside a subroutine. */

  if ((q = mxGetPr(OP)) == NULL)
    mexErrMsgTxt("mxGetPr(OP) failed");

  /* If both MU_A and MU_S specified, the matrix is twice as large */

  if (q[MU_A] != 0 && q[MU_S] != 0) 
    nAVox = 2*nVox;
  else
    nAVox =   nVox;

  plhs[0] = mxCreateDoubleMatrix(nMeas, nAVox, mxREAL);
  Atmpm   = mxCreateDoubleMatrix(nAVox, 1, mxREAL);

  A    = mxGetPr(plhs[0]);
  Atmp = mxGetPr(Atmpm);

  /* Initialize memory - zeros everywhere. */

  for (i = 0; i < nMeas*nAVox; i++)
    A[i] = 0.0;

  for (i = 0; i < nAVox; i++)
    Atmp[i] = 0.0;

  /* Do the actual computation */

  if (debug)
    mexPrintf("Meas: \n");

  for (iMeas = 0; iMeas < nMeas; iMeas++)
    {
      int j;

      TD3pt(SD, Medium, MeasList, OP, rvox, nVox, iMeas, Atmp, debug);

      if (debug)
        mexPrintf(" %d\n", iMeas+1); /* there is no mexFlush() */

      mxAssert(nAVox == nVox || nAVox == 2*nVox, "Length mismatch");

      /* Copy back from temporary memory to final storage and rescale
       * matrix elements by the volume of each voxel.  Do the
       * comparison (nAVox==nVox) outside the loop for efficiency.
       */

      for (j = 0; j < nVox; j++)
        A[j * nMeas + iMeas] = Atmp[j] * vvox[j];

      if (nAVox == 2*nVox)
        for (j = nVox; j < 2*nVox; j++)
          A[j * nMeas + iMeas] = Atmp[j] * vvox[j-nVox];
    }

  /* Done with voxel vectors */

  mxFree(rvox);
  mxFree(vvox);

  mxDestroyArray(Atmpm);

  if (debug)
    mexPrintf("\n");

  return;
}

/* ******************************************************************** */

static void TD3pt(const mxArray *SD, const mxArray *Medium, 
                  const mxArray *MeasList, const mxArray *OptProp,
                  const double *rv, int nvox, int iMeas, double *A, int debug)
{
  struct SimParam sp;

  double sa, da, Tgate, Tdelay, dtsrc, dtdet, *q;
  int    iSrc, iDet, iWvl, iGate, iDelay, i;
  int    nSrc, nDet, nWvl, nMeas, navox;
  mxArray *p;

  if (!mxIsStruct(SD))
    mexErrMsgTxt("TD3pt illegal SD struct");

  if (!mxIsStruct(Medium))
    mexErrMsgTxt("TD3pt illegal Medium struct");

  if ((q = mxGetPr(OptProp)) == NULL)
    mexErrMsgTxt("mxGetPr(OptProp) failed");

  sp.muflag[MU_A] = (q[MU_A] != 0);
  sp.muflag[MU_S] = (q[MU_S] != 0);

  /* **************************************************************** */

  /* Unpack the Matlab data structures -- Indices into arrays */

  nMeas = mxGetM(MeasList);

  if ((q = mxGetPr(MeasList)) == NULL)
    mexErrMsgTxt("mxGetPr(MeasList) failed");

  iSrc   = (int)q[nMeas * OFFSET_SRC   + iMeas] - 1;
  iDet   = (int)q[nMeas * OFFSET_DET   + iMeas] - 1;
  iWvl   = (int)q[nMeas * OFFSET_WVL   + iMeas] - 1;
  iGate  = (int)q[nMeas * OFFSET_GATE  + iMeas] - 1;
  iDelay = (int)q[nMeas * OFFSET_DELAY + iMeas] - 1;

  /* Unpack the Matlab data structures -- individual fields */

  /* All I care about is the number of wavelengths */

  if ((q = getVector(SD, "Lambda", &nWvl)) == NULL)
    mexErrMsgTxt("getVector() returns NULL");

  /* Get optode positions and amplitudes */

  if ((q = getVector(SD, "SrcPos", &nSrc)) == NULL)
    mexErrMsgTxt("getVector() returns NULL");

  nSrc /= 3;		  /* Size of the columns */

  sp.rsrc[0] = q[iSrc + 0 * nSrc];
  sp.rsrc[1] = q[iSrc + 1 * nSrc];
  sp.rsrc[2] = q[iSrc + 2 * nSrc];

  sa = getField(SD, "SrcAmp", iSrc);

  if ((q = getVector(SD, "DetPos", &nDet)) == NULL)
    mexErrMsgTxt("getVector() returns NULL");

  nDet /= 3;		  /* Size of the columns */

  sp.rdet[0] = q[iDet + 0 * nDet];
  sp.rdet[1] = q[iDet + 1 * nDet];
  sp.rdet[2] = q[iDet + 2 * nDet];

  da = getField(SD, "DetAmp", iDet);

  /* Detector characteristics */

  Tgate   = getField(SD,    "TimeGateWidth", iGate);
  Tdelay  = getField(SD,    "TimeDelay",     iDelay);

  /* Optical properties */

  sp.musp = getField(Medium, "Muspo",        iWvl);
  sp.mua  = getField(Medium, "Muao",         iWvl);
  sp.n    = getField(Medium, "idxRefr",      iWvl);
  sp.v    = C_VACUUM / sp.n;

  /* Geometry */

  if ((p = mxGetField(Medium, 0, "Slab_Thickness")) != NULL)
    sp.slabz = getField(Medium, "Slab_Thickness", 0);
  else
    sp.slabz = 0.0;                /* Not always relevant */

  /* Compute detector start and stop times (source is implicitly T=0) */

  dtsrc = dtdet = 0.0;

  if ((p = mxGetField(SD, 0, "SrcOffset")) != NULL)
    if ((q = mxGetPr(p)) != NULL)
      {
        int nSOffset = mxGetM(p);

        if ((mxGetM(p) != nSrc) || (mxGetN(p) != nWvl))
          mexErrMsgTxt("SD.SrcOffset size mismatch\n");

	dtsrc = q[iSrc + iWvl*nSOffset];
      }

  if ((p = mxGetField(SD, 0, "DetOffset")) != NULL)
    if ((q = mxGetPr(p)) != NULL)
      {
        int nDOffset = mxGetM(p);

        if ((mxGetM(p) != nDet) || (mxGetN(p) != nWvl))
          mexErrMsgTxt("SD.DetOffset size mismatch\n");

	dtdet = q[iDet + iWvl*nDOffset];
      }

  /* Shift times so that the light pulse always arrives at T=0 */

  sp.T1 = Tdelay + dtdet - dtsrc;
  sp.T2 = sp.T1 + Tgate;

  sp.debug = debug;

  /* Compute each element in the vector A */

  navox = (sp.muflag[MU_A] && sp.muflag[MU_S]) ? 2*nvox : nvox;

  TD3GateInt(&sp, A, navox, rv, nvox);

  /* Rescale computed results by source and detector amplitudes */

  for (i = 0; i < navox; i++)
    A[i] *= sa * da;

  return;
}

