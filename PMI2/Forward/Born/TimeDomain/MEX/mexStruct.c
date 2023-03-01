/* Code to simplify working with Matlab structures */

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

#include <stdio.h>

#include "mex.h"

/* Extract a single field out of a Matlab structure.  Offset is the
 *   offset into the array of doubles to use (use offset = 0 for scalars) */

double getField(const mxArray *mstruct, const char *name, unsigned int offset)
{
  mxArray *p;
  double *q;

  mxAssert(mstruct != NULL, "getField(): mstruct is NULL");
  mxAssert(name    != NULL, "getField(): name is NULL");
  mxAssert(offset  >= 0,    "getField(): illegal offset");

  if (!mxIsStruct(mstruct))
    mexErrMsgTxt("Illegal mstructure pointer");

  if ((p = mxGetField(mstruct, 0, name)) == NULL)
    {
      char buffer[256];

      snprintf(buffer, 255, "%s - field not found", name);
      mexErrMsgTxt(buffer);
    }

  if ((q = mxGetPr(p)) == NULL)
    mexErrMsgTxt("mxGetPr(p) failed");

  return q[offset];
}

/* Extract a data vector out of a Matlab structure.  Sz is the length
 *  of the data vector returned (use sz=NULL if you don't care) */

double *getVector(const mxArray * mstruct, const char *name, int *sz)
{
  mxArray *p;
  double *q;

  mxAssert(mstruct != NULL, "getVector(): mstruct is NULL");
  mxAssert(name    != NULL, "getVector(): name is NULL");

  if (!mxIsStruct(mstruct))
    mexErrMsgTxt("Illegal mstructure pointer");

  if ((p = mxGetField(mstruct, 0, name)) == NULL)
    {
      char buffer[256];

      snprintf(buffer, 256, "%s - field not found", name);
      mexErrMsgTxt(buffer);
    }

  if ((q = mxGetPr(p)) == NULL)
    mexErrMsgTxt("mxGetPr(p) failed");

  if (sz != NULL)
    *sz = mxGetNumberOfElements(p);

  return q;
}
