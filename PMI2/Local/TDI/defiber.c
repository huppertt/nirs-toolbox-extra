/*
 * Given a point "inside" the fiber (x1,x2) and a threshold level, do
 * a flood fill to find all points in the simply-connected domain A of
 * points "in" the fiber and greater than the threshold and return
 * their average.  
 *
 * meas = defiber(data, x1, x2, threshold, debug);
 */

#include "mex.h"
#include <string.h>

static unsigned char *_mask = NULL;

static double floodfill(int *npt, double *data, int n1, int n2, 
                        int x1, int x2, double thresh, int debug);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *data, threshold, *p, Y;
  int     x1, x2, n1, n2, N, debug;

  if ((nrhs < 4) || (nrhs > 5))
    mexErrMsgTxt("Incorrect input paramter count");
  else if (nlhs > 2)
    mexErrMsgTxt("Too many output arguments");
  else if (nlhs < 1)
    mexErrMsgTxt("Too few output arguments");

  /* data = array of doubles */

  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("data must be an array of doubles");

  n1 = mxGetM(prhs[0]);
  n2 = mxGetN(prhs[0]);

  if ((p = mxGetPr(prhs[0])) == NULL)
    mexErrMsgTxt("Can't get data out of mxArray");

  data = p;

  /* x1 = integer coordinate */

  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("x1 must be an integer");

  if ((p = mxGetPr(prhs[1])) == NULL)
    mexErrMsgTxt("Can't get x1 out of mxArray");

  x1 = (int)p[0];

  /* x2 = integer coordinate */

  if (!mxIsDouble(prhs[2]))
    mexErrMsgTxt("x2 must be an integer");

  if ((p = mxGetPr(prhs[2])) == NULL)
    mexErrMsgTxt("Can't get x2 out of mxArray");

  x2 = (int)p[0];

  /* threshold = double */

  if (!mxIsDouble(prhs[3]))
    mexErrMsgTxt("threshold must be a double");

  if ((p = mxGetPr(prhs[3])) == NULL)
    mexErrMsgTxt("Can't get threshod out of mxArray");

  threshold = (int)p[0];

  /* debug = optional integer */

  debug = 0;

  if (nrhs > 4)
    {
      if (!mxIsDouble(prhs[4]))
        mexErrMsgTxt("debug must be an integer");

      if ((p = mxGetPr(prhs[4])) == NULL)
        mexErrMsgTxt("Can't get debug out of mxArray");

      debug = (int)p[0];
    }

  /* Allocate space for the mask and set to zero */

  if ((_mask = mxMalloc(n1 * n2 * sizeof(unsigned char))) == NULL)
    mexErrMsgTxt("Error allocating space for mask");
  else
    memset(_mask, 0, n1 * n2 * sizeof(unsigned char));

  /* Call the real fill function */

  if (data[x2*n1 + x1] < threshold)
    {
      if (debug)
        mexPrintf("Trivial case, returning 0\n");

      Y = 0.0;                              /* Trivial case */
    }
  else
    {
      if (debug)
        mexPrintf("floodfill(&N, data, %d, %d, %d, %d, %d, %d\n",
                  n1, n2, x1, x2, threshold, debug);

      Y = floodfill(&N, data, n1, n2, x1, x2, threshold, debug);
    }

  /* Pass back as LHS */

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  p = mxGetPr(plhs[0]);

  p[0] = Y;

  if (nlhs > 1)
    {
      int i, j;

      plhs[1] = mxCreateDoubleMatrix(n1, n2, mxREAL);
      p = mxGetPr(plhs[1]);

      for (i = 0; i < n2; i++)
        for (j = 0; j < n1; j++)
          {
            int k = i * n1 + j;

            p[k] = (double)_mask[k];
          }
    }
      
  mxFree(_mask);
  _mask = NULL;

  return;
}

/***********************************************************************
 *
 * The actuall flood-fill code, recursively sums up all points above
 * threshold, marking off the mask as it goes.
 *
 **********************************************************************/

static double floodfill(int *npt, double *data, int n1, int n2, 
                        int x1, int x2, double thresh, int debug)
{
  double y1, y2, y3, y4, y;
  int    m1, m2, m3, m4, x = x2*n1 + x1;

  if (x >= n1*n2)
    {
      if (debug)
        mexPrintf("-- passed end of table [%d,%d]\n", x1, x2);

      *npt = 0;
      return 0.0;
    }

  if (debug)
    mexPrintf("Filling from data[%d,%d]\n", x1, x2);

  /* It shouldn't be possible to go out of bounds or be null */

  mxAssert(npt != NULL && data != NULL, "NULL pointer");
  mxAssert((x1 >= 0) && (x1 < n1), "x1 out of range");
  mxAssert((x2 >= 0) && (x2 < n2), "x2 out of range");
  mxAssert(_mask != NULL, "Data mask is NULL");

  *npt = 0;

  if (_mask[x] > 0)
    {
      if (debug)
        mexPrintf("-- data[%d,%d] masked out\n", x1, x2);

      return 0.0;
    }

  if (data[x] < thresh)
    {
      if (debug)
        mexPrintf("-- data[%d,%d] below threshold\n", x1, x2);

      return 0.0;
    }

  /* This is a valid point, add it in and mask it off */

  *npt += 1;
  y = data[x];
  _mask[x] = 1;

  m1 = m2 = m3 = m4 = 0;
  y1 = y2 = y3 = y4 = 0;

  if (x1 >= 1)
    y1 = floodfill(&m1, data, n1, n2, x1 - 1, x2, thresh, debug);

  if (x1 < n1-1)
    y2 = floodfill(&m2, data, n1, n2, x1 + 1, x2, thresh, debug);

  if (x2 >= 1)
    y3 = floodfill(&m3, data, n1, n2, x1, x2 - 1, thresh, debug);

  if (x2 < n2 - 1)
    y4 = floodfill(&m4, data, n1, n2, x1, x2 + 1, thresh, debug);

  y   += y1 + y2 + y3 + y4;
  npt += m1 + m2 + m3 + m4;

  return y;
}

