/*
 * Turn Medium.CompVol into a list of voxel centers and their
 * individual volumes.
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

#include "mex.h"

#include <math.h>
#include <string.h>
#include "mex.h"

/* Fight numerical round-off at the 1e-16 level */
#define FUDGE 1e-12

static int getUniform (const mxArray *, double **, double **);
static int getComputed(const mxArray *, double **, double **);
static int getList    (const mxArray *, double **, double **);

extern double *getVector(const mxArray *, const char *, int *);
extern double  getField(const mxArray *,  const char *, unsigned int);

/*
 * p must point to the contents of Medium.CompVol at start.
 * rv and vv are output parameters
 */

int getSampleVolume(const mxArray *p, double **rv, double **vv)
{
  const mxArray *q;
  int nvox = -1;

  *rv = NULL;
  *vv = NULL;

  if ((q = mxGetField(p, 0, "Type")) == NULL)
    mexErrMsgTxt("Medium.CompVol.Type not found");
  else
    {
      char buffer[1024];
      int n = 1024;

      if (mxGetString(q, buffer, n) == 1)
        mexErrMsgTxt("Medium.CompVol.Type not a string");

      /* Based on Medium.CompVol.Type, call the approprate subroutine */

#if (DEBUG)
      mexPrintf("Type string is \"%s\"\n", buffer);
#endif

      if      (strcasecmp(buffer, "uniform")  == 0)
        nvox = getUniform(p, rv, vv);
      else if (strcasecmp(buffer, "computed") == 0)
        nvox = getComputed(p, rv, vv);
      else if (strcasecmp(buffer, "list")     == 0)
        nvox = getList(p, rv, vv);
      else
        {
          mexPrintf("Unknown value for Medium.CompVol.Type, '%s'\n", buffer);
          mexErrMsgTxt("Unable to proceed");
        }
    }

  return nvox;
}

static int getUniform (const mxArray *p, double **rvox, double **vvox)
{
  double *X, *Y, *Z, *rv, *vv, sx, sy, sz;
  int    i, j, k, nx, ny, nz, nVox=-1;

  X = getVector(p, "X", &nx);
  Y = getVector(p, "Y", &ny);
  Z = getVector(p, "Z", &nz);
      
  nVox = nx * ny * nz;

  if (nVox == 0)
    mexErrMsgTxt("Medium.CompVol - volume has zero size");

  if ((rv = mxMalloc(3*nx*ny*nz*sizeof(double))) == NULL)
    mexErrMsgTxt("Malloc() of position vector failed");

  /* Generate table of (vector) position of every voxel, needed to
   * generate the forward matrx.
   *
   * In matlab, to which we need to maintain compatibility, Y (and
   * not X) is the fastest changing dimension.  
   */

  for (k = 0; k < nz; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        {
          int idx = j + ny*(i + nx*k);

          rv[3*idx + 0] = X[i];
          rv[3*idx + 1] = Y[j];
          rv[3*idx + 2] = Z[k];
        } 

  sx = getField(p, "XStep", 0);
  sy = getField(p, "YStep", 0);
  sz = getField(p, "ZStep", 0);

  if ((vv = mxMalloc(nVox*sizeof(double))) == NULL)
    mexErrMsgTxt("Malloc() of volume vector failed");

  for (i = 0; i < nVox; i++)
    vv[i] = sx*sy*sz;

  /* Done, set pointers for return */

  *rvox = rv;
  *vvox = vv;

  return nVox;
}

static int getComputed(const mxArray *p, double **rvox, double **vvox)
{
  const mxArray *q, *r;
  double x1, x2, y1, y2, z1, z2;
  double sx, sy, sz, *rv, *vv;
  int    i, j, k, nx, ny, nz, nVox;

  sx = getField(p, "XStep", 0);
  sy = getField(p, "YStep", 0);
  sz = getField(p, "ZStep", 0);

  if ((sx <= 0) || (sy <= 0) || (sz <= 0))
    mexErrMsgTxt("Step sizes must be strictly positive");

  /* The ranges can be specified as Medium.CompVol.X(1:2) or
   *  as Medium.CompVol.X1,Medium.CompVol.X2, etc.
   */

  if (((q = mxGetField(p, 0, "X1")) != NULL) &&
      ((r = mxGetField(p, 0, "X2")) != NULL))
    {
      /* Now that I know they exist... */
      double *R1 = getVector(p, "X1", NULL);
      double *R2 = getVector(p, "X2", NULL);

      x1 = R1[0];
      x2 = R2[0];
    }
  else
    {
      /* If this doesn't exist, it's a genuine error */
      double *R = getVector(p, "X", &nx);

      if (R == NULL)
        mexErrMsgTxt("Medium.CompVol.X not found");

      if (nx != 2)
        mexErrMsgTxt("Medium.CompVol.X does not have length of 2");
      
      x1 = R[0];
      x2 = R[1];
    }

  for (nx = 0; x2 >= x1 + nx*sx - FUDGE; nx++)  /* Count up voxels */
    ;

  /* Y */

  if (((q = mxGetField(p, 0, "Y1")) != NULL) &&
      ((r = mxGetField(p, 0, "Y2")) != NULL))
    {
      /* Now that I know they exist... */
      double *R1 = getVector(p, "Y1", NULL);
      double *R2 = getVector(p, "Y2", NULL);

      y1 = R1[0];
      y2 = R2[0];
    }
  else
    {
      /* If this doesn't exist, it's a genuine error */
      double *R = getVector(p, "Y", &ny);

      if (R == NULL)
        mexErrMsgTxt("Medium.CompVol.Y not found");

      if (ny != 2)
        mexErrMsgTxt("Medium.CompVol.Y does not have length of 2");
      
      y1 = R[0];
      y2 = R[1];
    }

  for (ny = 0; y2 >= y1 + ny*sy - FUDGE; ny++)  /* Count up voxels */
    ;

  /* Z */

  if (((q = mxGetField(p, 0, "Z1")) != NULL) &&
      ((r = mxGetField(p, 0, "Z2")) != NULL))
    {
      /* Now that I know they exist... */
      double *R1 = getVector(p, "Z1", NULL);
      double *R2 = getVector(p, "Z2", NULL);

      z1 = R1[0];
      z2 = R2[0];
    }
  else
    {
      /* If this doesn't exist, it's a genuine error */
      double *R = getVector(p, "Z", &nz);

      if (R == NULL)
        mexErrMsgTxt("Medium.CompVol.Z not found");

      if (nz != 2)
        mexErrMsgTxt("Medium.CompVol.Z does not have length of 2");
      
      z1 = R[0];
      z2 = R[1];
    }

  for (nz = 0; z2 >= z1 + nz*sz - FUDGE; nz++)  /* Count up voxels */
    ;

  /* Build the position vector */

  nVox = nx * ny * nz;

  if (nVox == 0)
    mexErrMsgTxt("Medium.CompVol - volume has zero size");

  if ((rv = mxMalloc(3*nx*ny*nz*sizeof(double))) == NULL)
    mexErrMsgTxt("Malloc() of position vector failed");

  /* Generate table of (vector) position of every voxel, needed to
   * generate the forward matrx.
   *
   * In matlab, to which we need to maintain compatibility, Y (and
   * not X) is the fastest changing dimension.  
   */

  for (k = 0; k < nz; k++)
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        {
          int idx = j + ny*(i + nx*k);

          rv[3*idx + 0] = x1 + i * sx;
          rv[3*idx + 1] = y1 + j * sy;
          rv[3*idx + 2] = z1 + k * sz;
        } 

  /* Build up the voxel size vector */

  if ((vv = mxMalloc(nx*ny*nz*sizeof(double))) == NULL)
    mexErrMsgTxt("Malloc() of volume vector failed");

  for (i = 0; i < nx*ny*nz; i++)
    vv[i] = sx*sy*sz;

  /* Done, set pointers for return */

  *rvox = rv;
  *vvox = vv;

  return nVox;
}

static int getList(const mxArray *p, double **rvox, double **vvox)
{
  double *X, *Y, *Z, *V, *rv, *vv;
  int     i, nx, ny, nz, nv;

  if ((X = getVector(p, "X", &nx)) == NULL)
    mexErrMsgTxt("Medium.CompVol.X not defined");
  if ((Y = getVector(p, "Y", &ny)) == NULL)
    mexErrMsgTxt("Medium.CompVol.Y not defined");
  if ((Z = getVector(p, "Z", &nz)) == NULL)
    mexErrMsgTxt("Medium.CompVol.Z not defined");
  if ((V = getVector(p, "Volume", &nv)) == NULL)
    mexErrMsgTxt("Medium.CompVol.Volume not defined");

  if ((nx != ny) || (nx != nz) || (nx != nv))
    mexErrMsgTxt("Medium.CompVol.* vector lengths must be the same");

  if ((rv = mxMalloc(3*nx*sizeof(double))) == NULL)
    mexErrMsgTxt("Malloc() of position vector failed");

  if ((vv = mxMalloc(  nx*sizeof(double))) == NULL)
    mexErrMsgTxt("Malloc() of position vector failed");

  /* Copy out of Medium.CompVol into my workspace */

  for (i = 0; i < nv; i++)
    {
      rv[3*i + 0] = X[i];
      rv[3*i + 1] = Y[i];
      rv[3*i + 2] = Z[i];

      vv[i]       = V[i];
    }

  /* Done, set pointers for return */

  *rvox = rv;
  *vvox = vv;

  return nv;
}
