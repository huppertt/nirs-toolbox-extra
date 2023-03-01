/* ---------------------------------------------------------------------
 *
 * Compute 3-point Green's Function, G(rv-rs,t'-0) * G(rd-rv, t-t'), for
 *  scattering and absorbing perturbations.  Semi-infinite medium.
 *
 * ------------------------------------------------------------------ */

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

extern double   tdgf(const double *, const double *,
                     double, const struct SimParam *);
extern double *tdggf(const double *, const double *,
                     double, const struct SimParam *);

extern struct SimParam *sp;

static double  sharedint(double, double, int);
static double  slab_phiA(double[], double[], double);
static double *slab_phiS(double[], double[], double);

static void   MoveSrc(const struct SimParam *, double *, double *);

extern double tdAAnalytBase(double dt, double rs[3], double rd[3]);
extern double tdSAnalytBase(double dt, double rs[3], double rd[3]);

/* --- Numerical integration routines ------------------------------- */

double tdAintegrand(double tvox, void *p)
{
  double tdet = *(double *)p;

  if (tvox < 0)
    return 0;
  else
    return sharedint(tvox, tdet, DOMUA);
}

double tdSintegrand(double tvox, void *p)
{
  double tdet = *(double *)p;

  if (tvox < 0)
    return 0;
  else
    return sharedint(tvox, tdet, DOMUS);
}

/* ------------------------------------------------------------------ */

static double sharedint(double tmid, double tdet, int f)
{
  double rsrc[3], rvox[3], rdet[3], phi = -1;

  /* Make local copies */

  memcpy(rsrc, sp->rsrc, 3*sizeof(double));
  memcpy(rdet, sp->rdet, 3*sizeof(double));
  memcpy(rvox, sp->rvox, 3*sizeof(double));

  /* 
   * Move sources/detectors 1 mean free path into medium.  Implicitly
   * assumes sources start on an interface.
   */

  MoveSrc(sp, rsrc, rdet);

  /* Calculate the fluence (abs) or the gradient of the fluence (scat) */

  if (f == DOMUA)
    {
      double phisrc, phidet;

      phisrc = slab_phiA(rsrc, rvox, tmid - 0.00);
      phidet = slab_phiA(rdet, rvox, tdet - tmid);

      /* Turn photon densities into fluences */

      phi = (sp->v * phisrc) * (sp->v * phidet);
    }
  else if (f == DOMUS)
    {
      double phisrc[3], phidet[3], *phip;

      /* Sum the image charges _then_ take the dot product */

      phip = slab_phiS(rsrc, rvox, tmid - 0.00);
      mxAssert(phip != NULL, "slab_phiS returns NULL");
      memcpy(phisrc, phip, 3*sizeof(double));

      phip = slab_phiS(rdet, rvox, tdet - tmid);
      mxAssert(phip != NULL, "slab_phiS returns NULL");
      memcpy(phidet, phip, 3*sizeof(double));

      phi = DOT(phisrc, phidet);

      /* Turn photon densities into fluences */

      phi = sp->v * phi * sp->v;
    }
  else
    mexErrMsgTxt("Unknown flag f");

  return phi;
}

/* ----------------------------------------------------------------------
 * Compute fluence from a "source" to a given voxel after an interval
 *   dt.  Handle slab boundaries using method of images.
 * ------------------------------------------------------------------- */

static double slab_phiA(double rs[3], double rv[3], double dt)
{
  double z1, r[3], zi;
  double phi, phi_src, phi_img;

  phi_src = phi_img = 0;

  memcpy(r, rs, 3*sizeof(double));

  /* One extrapolated boundary.  The physical boundary is implcitly
     assumed to be located at z=0! */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");

  if (sp->slabz > 0)
    z1 = -sp->zext;
  else
    z1 =  sp->zext;

  /* Field due to source charge */

  phi_src = tdgf(r, rv, dt, sp);

  if (phi_src <= 0)
    {
      /* This can happen if two points are far away and dt is small.
       * If phi_src is zero, no number of image charges will ever make
       * it positive.
       */

      return 0.0;
    }

  /* Only have one image charge to worry about */

  zi = GET_IMAGE(r[2], z1);

  /* Calculate fluence due to image charges */

  r[2] = zi;
  phi_img = -1 * tdgf(r, rv, dt, sp);

  /* Sum source and image to get fluence */

  phi = phi_src + phi_img;

  return phi;
}

/* ----------------------------------------------------------------------
 * Compute fluence from a "source" to a given voxel after an interval
 *   dt.  Handle slab boundaries using method of images.
 * ------------------------------------------------------------------- */

static double *slab_phiS(double rs[3], double rv[3], double dt)
{
  static double phi[3];
  double phi_src[3], phi_img[3], *pphi;
  double z1, r[3], zi;

  phi_src[0] = phi_src[1] = phi_src[2] = 0;
  phi_img[0] = phi_img[1] = phi_img[2] = 0;

  memcpy(r, rs, 3*sizeof(double));

  /* One extrapolated boundaries.  The physical boundary is implcitly
     assumed to be located at z=0! */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");

  if (sp->slabz > 0)
    z1 = -sp->zext;
  else
    z1 =  sp->zext;

  /* Field due to source charge */

  pphi = tdggf(r, rv, dt, sp);
  memcpy(phi_src, pphi, 3*sizeof(double));

  if ((SQR(phi_src[0]) + SQR(phi_src[1]) + SQR(phi_src[2])) == 0)
    {
      /* This can happen if two points are far away and dt is small.
       * If phi_src is zero, no number of image charges will ever change that
       */

      phi[0] = phi[1] = phi[2] = 0.0;

      return phi;
    }

  /* Only have one image charge to worry about */

  zi = GET_IMAGE(r[2], z1);

  /* Calculate fluence due to image charges */

  r[2] = zi;
  pphi = tdggf(r, rv, dt, sp);

  phi_img[0] = -1 * pphi[0];
  phi_img[1] = -1 * pphi[1];
  phi_img[2] = -1 * pphi[2];

  /* Sum source and image to get fluence, copy to static memory and
     pass back address */

  phi[0] = phi_src[0] + phi_img[0];
  phi[1] = phi_src[1] + phi_img[1];
  phi[2] = phi_src[2] + phi_img[2];
  
  return phi;
}

/* --- Analytical expressions --------------------------------------- */

double tdAanalytic(double dt, void *p)
{
  double rs[3], rsi[3], rd[3], rdi[3];
  double z1=0, phi1, phi2, phi3, phi4;

  /* Move source/detectors into medium */

  memcpy(rs, sp->rsrc, 3*sizeof(double));
  memcpy(rd, sp->rdet, 3*sizeof(double));

  MoveSrc(sp, rs, rd);

  /* Generate the position vectors for the image charges */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");
  else
    if (sp->slabz > 0)
      z1 = -sp->zext;
    else
      z1 =  sp->zext;

  memcpy(rsi, rs, 3*sizeof(double));  rsi[2] = GET_IMAGE(rs[2], z1);
  memcpy(rdi, rd, 3*sizeof(double));  rdi[2] = GET_IMAGE(rd[2], z1);

  /* phi1 -> neither are images.  Multiply by extra sp->v to make this
   *  agree with integrated form.  Will be divided out later. */

  phi1 = tdAAnalytBase(dt, rs , rd );

  /* phi2 -> source is image */

  phi2 = tdAAnalytBase(dt, rsi, rd );

  /* phi3 -> detector is image */

  phi3 = tdAAnalytBase(dt, rs , rdi);

  /* phi4 -> both are images */

  phi4 = tdAAnalytBase(dt, rsi, rdi);

  /* Semi-infinite is simple, sum charges directly */

  return(phi1 - phi2 - phi3 + phi4);
}

double tdSanalytic(double tdet, void *p)
{
  double rs[3], rsi[3], rd[3], rdi[3];
  double z1=0, phi1, phi2, phi3, phi4;

  /* Move source/detectors into medium */

  memcpy(rs, sp->rsrc, 3*sizeof(double));
  memcpy(rd, sp->rdet, 3*sizeof(double));

  MoveSrc(sp, rs, rd);

  /* Generate the position vectors for the image charges */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");
  else
    if (sp->slabz > 0)
      z1 = -sp->zext;
    else
      z1 =  sp->zext;

  memcpy(rsi, rs, 3*sizeof(double));  rsi[2] = GET_IMAGE(rs[2], z1);
  memcpy(rdi, rd, 3*sizeof(double));  rdi[2] = GET_IMAGE(rd[2], z1);

  /* phi1 -> neither are images.  Multiply by extra sp->v to make this
   *  agree with integrated form.  Will be divided out later. */

  phi1 = tdSAnalytBase(tdet, rs, rd);

  /* phi2 -> source is image */

  phi2 = tdSAnalytBase(tdet, rsi, rd);

  /* phi3 -> detector is image */

  phi3 = tdSAnalytBase(tdet, rs, rdi);

  /* phi4 -> both are images */

  phi4 = tdSAnalytBase(tdet, rsi, rdi);

  /* Semi-infinite is simple */

  return(phi1 - phi2 - phi3 + phi4);
}

/* ****************************************************************** */

static void MoveSrc(const struct SimParam *sp, double *rs, double *rd)
{
  /* The boundary is implicitly at zero.  If slabz is defined, we use
   * its sign to decide which side of the boundary is zero, otherwise
   * we look at the detectors.  If the detectors are also at zero,
   * assume the positive direction is tissue. 
   */

  if      (sp->slabz < 0)
    {
      /* Move optodes to interface from outside volume */

      if (rs[2] >  0) rs[2] = 0;
      if (rd[2] >  0) rd[2] = 0;

      if (rs[2] == 0) rs[2] -= 1/sp->musp;
      if (rd[2] == 0) rd[2] -= 1/sp->musp;
    }
  else if (sp->slabz > 0)
    {
      if (rs[2] <  0) rs[2] = 0;
      if (rd[2] <  0) rd[2] = 0;

      if (rs[2] == 0) rs[2] += 1/sp->musp;
      if (rd[2] == 0) rd[2] += 1/sp->musp;
    }
  else
    mexErrMsgTxt("Illegal zero-thickness slab");

  return;
}

