/* ***********************************************************************
 *
 * Compute two-point Green's functions (Time-Domain DOT)
 *
 * From Matlab:
 *   Phi0 = mextd2ptXX(SD, Medium, [MeasList], [Debug]);
 *
 * Jonathan Stott, 2002
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

/* ******************************************************************** */

static double TD2pt(const mxArray *, 
                    const mxArray *, const mxArray *, int, int);

extern double *getVector(const mxArray *, const char *, int *);
extern double   getField(const mxArray *,  const char *, unsigned int);
extern double TD2GateInt(struct SimParam *);

/* ******************************************************************** */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *SD = NULL, *Medium = NULL, *MeasList = NULL;
  int i, nMeas = 0, debug;
  double *phi;

  if ((nrhs < 2) || (nrhs > 4))
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

  /* Look for the debug flag */

  debug = 0;

  if (nrhs == 4)
    {
      if (!mxIsDouble(prhs[3]))
        mexErrMsgTxt("Debug flag must be an integer");
      else
        {
          double *p = mxGetPr(prhs[3]);

          if (p != NULL)
            debug = (p[0] != 0.0);
        }
    }

  /* Create matrix for the return argument. */

  plhs[0] = mxCreateDoubleMatrix(nMeas, 1, mxREAL);

  phi = mxGetPr(plhs[0]);

  /* Initialize memory */

  for (i = 0; i < nMeas; i++)
    phi[i] = 0.0;

  /* Do the actual computation; loop over measurements */

  if (debug)
    mexPrintf("Meas: ");

  for (i = 0; i < nMeas; i++)
    {
      phi[i] = TD2pt(SD, Medium, MeasList, i, debug);

      if (debug)
        mexPrintf(" %d", i+1);
    }

  if (debug)
    mexPrintf("\n");

  return;
}

/* ******************************************************************** */

static double TD2pt(const mxArray *SD, const mxArray *Medium, 
                    const mxArray *MeasList, int iMeas, int debug)
{
  double sa, da, gatewidth, delaytime, phi, *q, dtsrc, dtdet;
  int    nMeas, iSrc, iDet, iWvl, iGate, iDelay, nSrc, nDet, nWvl;
  struct SimParam sp;
  mxArray *p;

  if (!mxIsStruct(SD))
    mexErrMsgTxt("TD2pt illegal SD struct");

  if (!mxIsStruct(Medium))
    mexErrMsgTxt("TD2pt illegal Medium struct");

  /* Unpack the Matlab data structures -- Indices into arrays */

  nMeas  = mxGetM(MeasList);

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

  nSrc /= 3;  /* size of first column */

  sp.rsrc[0] = q[nSrc * 0 + iSrc];
  sp.rsrc[1] = q[nSrc * 1 + iSrc];
  sp.rsrc[2] = q[nSrc * 2 + iSrc];

  sa = getField(SD, "SrcAmp", iSrc);

  if ((q = getVector(SD, "DetPos", &nDet)) == NULL)
    mexErrMsgTxt("getVector() returns NULL");

  nDet /= 3;  /* size of first column */

  sp.rdet[0] = q[nDet * 0 + iDet];
  sp.rdet[1] = q[nDet * 1 + iDet];
  sp.rdet[2] = q[nDet * 2 + iDet];

  da = getField(SD, "DetAmp", iDet);

  /* Temporal characteristics */

  gatewidth = getField(SD, "TimeGateWidth", iGate);
  delaytime = getField(SD, "TimeDelay",     iDelay);

  /* Optical properties */

  sp.musp= getField(Medium, "Muspo",   iWvl);
  sp.mua = getField(Medium, "Muao",    iWvl);
  sp.n   = getField(Medium, "idxRefr", iWvl);
  sp.v   = C_VACUUM / sp.n;

  /* Geometry */

  if ((p = mxGetField(Medium, 0, "Slab_Thickness")) != NULL)
    sp.slabz = getField(Medium, "Slab_Thickness", 0);
  else
    sp.slabz = 0.0;

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

  /* Shift times so the the light pulse always arrives at T=0 */
  /*  mexPrintf("dt = %f ns\n", (dtdet - dtsrc)*1e9); */
  sp.T1 = delaytime + dtdet - dtsrc;
  sp.T2 = sp.T1 + gatewidth;

  sp.debug = debug;

  /* 
   * Now that I've extracted the data out of the Matlab structures,
   * call the GSL integration routines.  Pass sp as a pointer just
   * to save memory/time.
   */

  phi = sa * da * TD2GateInt(&sp);

  return phi;
}
