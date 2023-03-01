/* ---------------------------------------------------------------------
 * quadpack routines return 0 on success, print error messag if nonzero
 * --------------------------------------------------------------------- */

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

#include "td2pt.h"
#include "gsl/gsl_errno.h"

void ierrmsg(int err, const char *func, double abserr);

/* Replace the default error handler */

void ierrhnd(const char *msg, const char *file, int line, int ierr)
{
  if (ierr == 0)
    {
      char buffer[1024];

      snprintf(buffer, 1023, 
               "Unknown error in %s line %d: %s", file, line, msg);
      mexWarnMsgTxt(buffer);
    }
  else
    {
#if 1
      mexPrintf("%s[%d]: %s\n", file, line, msg);
#else
      ierrmsg(ierr, file, 0.0);
#endif
    }

  return;
}

/* Print a nice error message */

void ierrmsg(int err, const char *func, double abserr)
{
  if (err != GSL_SUCCESS)
    {
      mexPrintf("%s: %s [abserr = %e]\n", func, gsl_strerror(err), abserr);

#if 1
      mexErrMsgTxt("GSL library error");
#endif
    }

  return;
}
