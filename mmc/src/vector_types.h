/***************************************************************************//**
**  \mainpage Mesh-based Monte Carlo (MMC) - a 3D photon simulator
**
**  \author Qianqian Fang <q.fang at neu.edu>
**  \copyright Qianqian Fang, 2010-2018
**
**  \section sref Reference:
**  \li \c (\b Fang2010) Qianqian Fang, <a href="http://www.opticsinfobase.org/abstract.cfm?uri=boe-1-1-165">
**          "Mesh-based Monte Carlo Method Using Fast Ray-Tracing 
**          in Plücker Coordinates,"</a> Biomed. Opt. Express, 1(1) 165-175 (2010).
**  \li \c (\b Fang2012) Qianqian Fang and David R. Kaeli, 
**           <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-3-12-3223">
**          "Accelerating mesh-based Monte Carlo method on modern CPU architectures,"</a> 
**          Biomed. Opt. Express 3(12), 3223-3230 (2012)
**  \li \c (\b Yao2016) Ruoyang Yao, Xavier Intes, and Qianqian Fang, 
**          <a href="https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-1-171">
**          "Generalized mesh-based Monte Carlo for wide-field illumination and detection 
**           via mesh retessellation,"</a> Biomed. Optics Express, 7(1), 171-184 (2016)
**
**  \section slicense License
**          GPL v3, see LICENSE.txt for details
*******************************************************************************/

/***************************************************************************//**
\file    vector_types.h

\brief   Definitions of the basic short vector data structures
*******************************************************************************/

#ifndef _MMC_VECTOR_H
#define _MMC_VECTOR_H

/**
 \struct MMC_float4 vector_types.h
 \brief  floating-point quadraplet {x,y,z,w}

 the data structure is 16byte aligned to facilitate SSE operations
*/

typedef struct MMC_float4{
    float x,y,z,w;
} float4 __attribute__ ((aligned(16)));

/**
 \struct MMC_float3 vector_types.h
 \brief  floating-point triplet {x,y,z}

 if SSE is enabled, float3 is identical to float4
*/

#ifdef MMC_USE_SSE
 typedef struct MMC_float4 float3;
#else
 typedef struct MMC_float3{
    float x,y,z;
 } float3;
#endif

/**
 \struct MMC_int2 vector_types.h
 \brief  integer pair {ix,iy}
*/

typedef struct MMC_int2{
    int x,y;
} int2;

/**
 \struct MMC_int3 vector_types.h
 \brief  integer triplet {ix,iy,iz}
*/

typedef struct MMC_int3{
    int x,y,z;
} int3;

/**
 \struct MMC_int4 vector_types.h
 \brief  unsigned integer quadraplet {ix,iy,iz,iw}
*/
typedef struct MMC_int4{
    int x,y,z,w;
} int4;

/**
 \struct MMC_uint3 vector_types.h
 \brief  unsigned integer triplet {ix,iy,iz}
*/
typedef struct MMC_uint3{
    unsigned int x,y,z;
} uint3;

/**
 \struct MMC_uint2 vector_types.h
 \brief  unsigned integer pair {ix,iy}
*/
typedef struct MMC_uint2{
    unsigned int x,y;
} uint2;

#endif
