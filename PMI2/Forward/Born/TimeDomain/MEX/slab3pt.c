/* ---------------------------------------------------------------------
 *
 * Compute 3-point Green's Function, G(rv-rs,t'-0) * G(rd-rv, t-t'), for
 *  scattering and absorbing perturbations.  Slab boundaries.
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

/* ------------------------------------------------------------------ */

double  tdAintegrand(double tvox, void *p)
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
  double z1, z2, chg, r[3], zi1, zi2, zii1, zii2;
  double phi_src, phi_img, phi_img1, phi_img2, dphi, phi;
  int    i;

  memcpy(r, rs, 3*sizeof(double));

  /* Two extrapolated boundaries.  One physical boundary is implcitly
     assumed to be located at z=0! */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");

  if (sp->slabz > 0)
    {
      z1 = 0         - sp->zext;
      z2 = sp->slabz + sp->zext;
    }
  else
    {
      z1 = sp->slabz - sp->zext;
      z2 = 0         + sp->zext;
    }

  /* Field due to source charge */

  phi_src = tdgf(r, rv, dt, sp);
  phi_img = 0;
  chg     = 1;

  if (phi_src <= 0)
    {
      /* This can happen if two points are far away and dt is small.
       * If phi_src is zero, no number of image charges will ever change that
       */

      return 0.0;
    }

  /* Initialize image depths to source position */

  zi1 = r[2];
  zi2 = r[2];

  /* Loop over image charges, maximum of MAXIMG images */

  for (i = 0; i < MAXIMG; i++)
    {
      zii1 = GET_IMAGE(zi2, z1);
      zii2 = GET_IMAGE(zi1, z2);

      /* Calculate fluence due to image charges */

      chg = -chg;

      r[2] = zii1;
      phi_img1 = chg * tdgf(r, rv, dt, sp);

      r[2] = zii2;
      phi_img2 = chg * tdgf(r, rv, dt, sp);

      /* Calculate relative change */

      dphi = ABS(phi_img1 + phi_img2) / ABS(phi_src + phi_img);

      /* Record this last iteration */

      phi_img += phi_img1 + phi_img2;

      /* Stop if we've converged sufficiently */

      if (dphi < PHI_TOL)
        break;
      else
        {
          /* Prepare for next iteration */

          zi1 = zii1;
          zi2 = zii2;
        }
    }

  if (i == MAXIMG)
    mexErrMsgTxt("Fluence failed to converge");
#if (DEBUG)
  else
    mexPrintf("Used %d image charges (A)\n", 2*(i+1));
#endif

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
  double phi_src[3], phi_img[3], phi_img1[3], phi_img2[3], *pphi;
  double z1, z2, chg, r[3], zi1, zi2, zii1, zii2, dphi;
  int    i;

  memcpy(r, rs, 3*sizeof(double));

  /* Two extrapolated boundaries.  One physical boundary is implcitly
     assumed to be located at z=0! */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");

  if (sp->slabz > 0)
    {
      z1 = 0         - sp->zext;
      z2 = sp->slabz + sp->zext;
    }
  else
    {
      z1 = sp->slabz - sp->zext;
      z2 = 0         + sp->zext;
    }

  /* Field due to source charge */

  pphi = tdggf(r, rv, dt, sp);
  memcpy(phi_src, pphi, 3*sizeof(double));

  phi_img[0] = phi_img[1] = phi_img[2] = 0.0;
  chg        = 1;

  if ((SQR(phi_src[0]) + SQR(phi_src[1]) + SQR(phi_src[2])) == 0)
    {
      /* This can happen if two points are far away and dt is small.
       * If phi_src is zero, no number of image charges will ever change that
       */

      phi[0] = phi[1] = phi[2] = 0.0;

      return phi;
    }

  /* Initialize image depths to source position */

  zi1 = r[2];
  zi2 = r[2];

  /* Loop over image charges, maximum of MAXIMG images */

  for (i = 0; i < MAXIMG; i++)
    {
      zii1 = GET_IMAGE(zi2, z1);
      zii2 = GET_IMAGE(zi1, z2);

      /* Calculate fluence due to image charges */

      chg = -chg;

      r[2] = zii1;
      pphi = tdggf(r, rv, dt, sp);

      phi_img1[0] = chg * pphi[0];
      phi_img1[1] = chg * pphi[1];
      phi_img1[2] = chg * pphi[2];

      r[2] = zii2;
      pphi = tdggf(r, rv, dt, sp);

      phi_img2[0] = chg * pphi[0];
      phi_img2[1] = chg * pphi[1];
      phi_img2[2] = chg * pphi[2];

      /* Calculate relative change */

      dphi = (SQR(phi_img1[0] + phi_img2[0]) +
              SQR(phi_img1[1] + phi_img2[1]) + SQR(phi_img1[2] + phi_img2[2]));
      dphi /=(SQR(phi_src[0] + phi_img[0]) +
              SQR(phi_src[1] + phi_img[1]) + SQR(phi_src[2] + phi_img[2]));
      dphi = sqrt(dphi);

      /* Record this last iteration */

      phi_img[0] += phi_img1[0] + phi_img2[0];
      phi_img[1] += phi_img1[1] + phi_img2[1];
      phi_img[2] += phi_img1[2] + phi_img2[2];

      /* Stop if we've converged sufficiently */

      if (dphi < PHI_TOL)
        break;
      else
        {
          /* Prepare for next iteration */
          
          zi1 = zii1;
          zi2 = zii2;
        }
    }

  if (i == MAXIMG)
    mexErrMsgTxt("Fluence failed to converge");
#if (DEBUG)
  else
    mexPrintf("Used %d image charges (S)\n", 2*(i+1));
#endif

  /* Copy to static memory and pass back address */
  phi[0] = phi_src[0] + phi_img[0];
  phi[1] = phi_src[1] + phi_img[1];
  phi[2] = phi_src[2] + phi_img[2];
  
  return phi;
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

      if (rs[2] > 0)         rs[2] = 0;
      if (rs[2] < sp->slabz) rs[2] = sp->slabz;

#if (MOVE_DETS)
      if (rd[2] > 0)         rd[2] = 0;
      if (rd[2] < sp->slabz) rd[2] = sp->slabz;
#endif

      if (rs[2] > sp->slabz/2)
        rs[2] -= 1/sp->musp;
      else
        rs[2] += 1/sp->musp;

#if (MOVE_DETS)
      if (rd[2] > sp->slabz/2)
        rd[2] -= 1/sp->musp;
      else
        rd[2] += 1/sp->musp;
#endif
    }
  else if (sp->slabz > 0)
    {
      if (rs[2] < 0)         rs[2] = 0;
      if (rs[2] > sp->slabz) rs[2] = sp->slabz;

#if (MOVE_DETS)
      if (rd[2] < 0)         rd[2] = 0;
      if (rd[2] > sp->slabz) rd[2] = sp->slabz;
#endif

      if (rs[2] < sp->slabz/2)
        rs[2] += 1/sp->musp;
      else
        rs[2] -= 1/sp->musp;

#if (MOVE_DETS)
      if (rd[2] < sp->slabz/2)
        rd[2] += 1/sp->musp;
      else
        rd[2] -= 1/sp->musp;
#endif
    }
  else
    mexErrMsgTxt("Illegal zero-thickness slab");

  return;
}

/* ****************************************************************** */

double tdAanalytic(double dt)
{
  double rs0[3], rd0[3], zs[MAXIMG][2], zd[MAXIMG][2], rsrc[3], rdet[3];
  double rsv, rvd, z1=0, z2=0, phi;
  int    i, j, nimg, charge;

  /* Generate the position vectors for the image charges */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");
  else
    if (sp->slabz > 0)
      {
        z1 = 0         - sp->zext;
        z2 = sp->slabz + sp->zext;
      }
    else
      {
        z1 = sp->slabz - sp->zext;
        z2 = 0         + sp->zext;
      }

  /* Move source/detectors into medium */

  memcpy(rs0, sp->rsrc, 3*sizeof(double));
  memcpy(rd0, sp->rdet, 3*sizeof(double));

  MoveSrc(sp, rs0, rd0);

  /* Initialize src/det vectors.  Compute image charge vectors only
     once.  Do all image charges now (GET_IMAGE() is not an expensive
     function call) */

  zs[0][0] = zs[0][1] = rs0[2];

  for (i = 1; i < MAXIMG; i++)
    {
      zs[i][0] = GET_IMAGE(zs[i-1][1], z1);
      zs[i][1] = GET_IMAGE(zs[i-1][0], z2);
    }

  zd[0][0] = zd[0][1] = rd0[2];

  for (j = 1; j < MAXIMG; j++)
    {
      zd[j][0] = GET_IMAGE(zd[j-1][1], z1);
      zd[j][1] = GET_IMAGE(zd[j-1][0], z2);
    }

  /* Initialize temporary vectors */

  memcpy(rsrc, rs0, 3*sizeof(double));
  memcpy(rdet, rd0, 3*sizeof(double));

  /* Initialize phi with the infinite-medium fluence */

  rsrc[2] = zs[0][0];
  rdet[2] = zd[0][0];
#if (DEBUG)
  mexPrintf("gen=%d, src_z=%f, det_z=%f, chg=%d\n", 0, rsrc[2], rdet[2], 1);
#endif

  rsv = SUBNORM(sp->rvox, rsrc);
  rvd = SUBNORM(rdet, sp->rvox);

  phi = tdAAnalytBase(dt, rsrc, rdet);

  /* For all orders of charges, i=0 and j=0 are the actual charges */

  for (nimg = 1; nimg < MAXIMG; nimg++)
    {
      double phiold = phi;

      /* This is not the most efficient way to count the additional
         elements in A, but it it at least demonstrably correct */

      for (i = -nimg; i <= nimg; i++)
        for (j = -nimg; j <= nimg; j++)
          {
            if ((ABS(i) != nimg) && (ABS(j) != nimg))
              continue;

            charge = (ABS(i+j) % 2) == 0 ? 1 : -1;

            /* Create appropriate position vectors */

            if (i >= 0)
              rsrc[2] = zs[ i][0];
            else
              rsrc[2] = zs[-i][1];

            if (j >= 0)
              rdet[2] = zd[ j][0];
            else
              rdet[2] = zd[-j][1];

#if (DEBUG)
            mexPrintf("gen=%d, src_z=%f, det_z=%f, chg=%d\n", 
                      nimg, rsrc[2], rdet[2], charge);
#endif

            rsv = SUBNORM(sp->rvox, rsrc);
            rvd = SUBNORM(rdet, sp->rvox);

            /* Fluence due to this distribution of "sources" */

            phi += charge * tdAAnalytBase(dt, rsrc, rdet);
          }

      /* Stop if we've converged sufficiently */

      if (ABS((phi-phiold)/(phi+phiold)) < PHI_TOL)
        break;

      /* If phi is still exactly zero after adding image charges, no
         number of additional image charges is going to alter that */
      
      if ((phi == 0.0) && (nimg > 0))
        {
#if (DEBUG)
          mexWarnMsgTxt("Phi below numeric precision, assuming zero");
#endif
          break;
        }
    }

  if (nimg == MAXIMG)
    {
      mexErrMsgTxt("Used all MAXIMG image charges and failed to converge\n");
      return 0.0;
    }
#if (DEBUG)
  else
    mexPrintf("Used %d image charges (A)\n", nimg);
#endif

  return phi;
}

double tdSanalytic(double dt)
{
  double rs0[3], rd0[3], zs[MAXIMG][2], zd[MAXIMG][2], rsrc[3], rdet[3];
  double rsv, rvd, z1=0, z2=0, phi;
  int    i, j, nimg, charge;

  /* Generate the position vectors for the image charges */

  if (sp->zext <= 0)
    mexErrMsgTxt("sp->zext is non-positive");
  else
    if (sp->slabz > 0)
      {
        z1 = 0         - sp->zext;
        z2 = sp->slabz + sp->zext;
      }
    else
      {
        z1 = sp->slabz - sp->zext;
        z2 = 0         + sp->zext;
      }

  /* Move source/detectors into medium */

  memcpy(rs0, sp->rsrc, 3*sizeof(double));
  memcpy(rd0, sp->rdet, 3*sizeof(double));

  MoveSrc(sp, rs0, rd0);

  /* Initialize src/det vectors.  Compute image charge vectors only
     once.  Do all image charges now (GET_IMAGE() is not an expensive
     function call) */

  zs[0][0] = zs[0][1] = rs0[2];

  for (i = 1; i < MAXIMG; i++)
    {
      zs[i][0] = GET_IMAGE(zs[i-1][1], z1);
      zs[i][1] = GET_IMAGE(zs[i-1][0], z2);
    }

  zd[0][0] = zd[0][1] = rd0[2];

  for (j = 1; j < MAXIMG; j++)
    {
      zd[j][0] = GET_IMAGE(zd[j-1][1], z1);
      zd[j][1] = GET_IMAGE(zd[j-1][0], z2);
    }

  /* Initialize temporary vectors */

  memcpy(rsrc, rs0, 3*sizeof(double));
  memcpy(rdet, rd0, 3*sizeof(double));

  /* Initialize phi with the infinite-medium fluence */

  rsrc[2] = zs[0][0];
  rdet[2] = zd[0][0];
#if (DEBUG)
  mexPrintf("gen=%d, src_z=%f, det_z=%f, chg=%d\n", 0, rsrc[2], rdet[2], 1);
#endif

  rsv = SUBNORM(sp->rvox, rsrc);
  rvd = SUBNORM(rdet, sp->rvox);

  phi = tdSAnalytBase(dt, rsrc, rdet);

  /* For all orders of charges, i=0 and j=0 are the actual charges */

  for (nimg = 1; nimg < MAXIMG; nimg++)
    {
      double phiold = phi;

      /* This is not the most efficient way to count the additional
         elements in A, but it it at least demonstrably correct */

      for (i = -nimg; i <= nimg; i++)
        for (j = -nimg; j <= nimg; j++)
          {
            if ((ABS(i) != nimg) && (ABS(j) != nimg))
              continue;

            charge = (ABS(i+j) % 2) == 0 ? 1 : -1;

            /* Create appropriate position vectors */

            if (i >= 0)
              rsrc[2] = zs[ i][0];
            else
              rsrc[2] = zs[-i][1];

            if (j >= 0)
              rdet[2] = zd[ j][0];
            else
              rdet[2] = zd[-j][1];

#if (DEBUG)
            mexPrintf("gen=%d, src_z=%f, det_z=%f, chg=%d\n", 
                      nimg, rsrc[2], rdet[2], charge);
#endif

            rsv = SUBNORM(sp->rvox, rsrc);
            rvd = SUBNORM(rdet, sp->rvox);

            /* Fluence due to this distribution of "sources" */

            phi += charge * tdSAnalytBase(dt, rsrc, rdet);
        }

      /* Stop if we've converged sufficiently */

      if (ABS((phi-phiold)/(phi+phiold)) < PHI_TOL)
        break;

      /* If phi is still exactly zero after adding image charges, no
         number of additional image charges is going to alter that */

      if ((phi == 0.0) && (nimg > 0))
        {
#if (DEBUG)
          mexWarnMsgTxt("Phi below numeric precision, assuming zero");
#endif
          break;
        }
    }

  if (nimg == MAXIMG)
    {
      mexErrMsgTxt("Used all MAXIMG image charges and failed to converge\n");
      return 0.0;
    }
#if (DEBUG)
  else
    mexPrintf("Used %d image charges (S)\n", nimg);
#endif

  return phi;
}
