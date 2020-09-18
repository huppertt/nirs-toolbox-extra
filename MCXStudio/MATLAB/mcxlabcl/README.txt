= MCXLAB-CL: MCX for MATLAB and GNU Octave =

Author: Qianqian Fang <q.fang at neu.edu>
License: GNU General Public License version 3 (GPLv3)
Version: this package is part of Monte Carlo eXtreme OpenCL (MCX-CL) v2018.3

<toc>


== # Introduction ==

MCXLAB-CL is the native MEX version of MCX for Matlab and GNU Octave. It compiles
the entire MCX code into a MEX function which can be called directly inside
Matlab or Octave. The input and output files in MCX are replaced by convenient
in-memory struct variables in MCXLAB-CL, thus, making it much easier to use 
and interact. Matlab/Octave also provides convenient plotting and data
analysis functions. With MCXLAB-CL, your analysis can be streamlined and speed-
up without involving disk files.

Because MCXLAB-CL contains the exact computational codes for the GPU calculations
as in the MCX binaries, MCXLAB-CL is expected to have identical performance when
running simulations. By default, we compile MCXLAB-CL with the support of recording
detected photon partial path-lengths (i.e. the "make det" option).

== # Installation ==

To download MCXLAB-CL, please visit [http://mcx.sourceforge.net/cgi-bin/index.cgi?Download#Download_the_Latest_Release this link]. 
If you choose to [http://mcx.sourceforge.net/cgi-bin/index.cgi?register/mcx register], 
you will have an option to be notified for any future updates.

The system requirements for MCXLAB-CL are the same as MCX: you have to make
sure that you have a CUDA-capable graphics card with properly configured 
CUDA driver (you can run the standard MCX binary first to test if your 
system is capable to run MCXLAB-CL). Of course, you need to have either Matlab
or Octave installed.

Once you set up the CUDA toolkit and NVIDIA driver, you can then add the 
"mcxlabcl" directory to your Matlab/Octave search path using the addpath command.
If you want to add this path permanently, please use the "pathtool" 
command, or edit your startup.m (~/.octaverc for Octave).

If everything works ok, typing "help mcxlabcl" in Matlab/Octave will print the
help information. If you see any error, particularly any missing libraries,
please make sure you have downloaded the matching version built for your
platform.


== # How to use MCXLAB-CL in MATLAB/Octave ==

To learn the basic usage of MCXLAB-CL, you can type

  help mcxlabcl

in Matlab/Octave to see the help information regarding how to use this 
function. The help information is listed below. You can find the input/output 
formats and examples. The input cfg structure has very similar field names as
the verbose command line options in MCX.

<pre> ====================================================================
       MCXLAB-CL - Monte Carlo eXtreme (MCX) for MATLAB/GNU Octave
 --------------------------------------------------------------------
        Copyright (c) 2018-2019 Qianqian Fang <q.fang at neu.edu>
                       URL: http://mcx.space
 ====================================================================
 
  Format:
     fluence=mcxlabcl(cfg);
        or
     [fluence,detphoton,vol]=mcxlabcl(cfg);
     [fluence,detphoton,vol]=mcxlabcl(cfg, option);
 
  Input:
     cfg: a struct, or struct array. Each element of cfg defines 
          the parameters associated with a simulation. 
          if cfg='gpuinfo': return the supported GPUs and their parameters,
          see sample script at the bottom
     option: (optional), options is a string, specifying additional options
          option='preview': this plots the domain configuration using mcxpreview(cfg)
          option='opencl':  force using mcxcl.mex* instead of mcx.mex* on NVIDIA/AMD/Intel hardware
          option='cuda':    force using mcx.mex* instead of mcxcl.mex* on NVIDIA GPUs
 
     if one defines USE_MCXCL=1 in MATLAB command line window, all following
     mcxlab and mcxlabcl calls will use mcxcl.mex; by setting option='cuda', one can
     force both mcxlab and mcxlabcl to use mcx (cuda version). Similarly, if
     USE_MCXCL=0, all mcxlabcl and mcxlab call will use mcx.mex by default, unless
     one set option='opencl'.
 
     cfg may contain the following fields:
 
 == Required ==
      *cfg.nphoton:    the total number of photons to be simulated (integer)
      *cfg.vol:        a 3D array specifying the media index in the domain
                       can be uint8, uint16, uint32, single or double
                       arrays.
      *cfg.prop:       an N by 4 array, each row specifies [mua, mus, g, n] in order.
                       the first row corresponds to medium type 0
                       (background) which is typically [0 0 1 1]. The
                       second row is type 1, and so on. The background
                       medium (type 0) has special meanings: a photon
                       terminates when moving from a non-zero to zero voxel.
      *cfg.tstart:     starting time of the simulation (in seconds)
      *cfg.tstep:      time-gate width of the simulation (in seconds)
      *cfg.tend:       ending time of the simulation (in second)
      *cfg.srcpos:     a 1 by 3 vector, the position of the source in grid unit
      *cfg.srcdir:     a 1 by 3 vector, specifying the incident vector; if srcdir
                       contains a 4th element, it specifies the focal length of
                       the source (only valid for focuable src, such as planar, disk,
                       fourier, gaussian, pattern, slit, etc); if the focal length
                       is nan, all photons will be launched isotropically regardless
 
 == MC simulation settings ==
       cfg.seed:       seed for the random number generator (integer) [0]
                       if set to a uint8 array, the binary data in each column is used 
                       to seed a photon (i.e. the "replay" mode)
       cfg.respin:     repeat simulation for the given time (integer) [1]
       cfg.isreflect:  [1]-consider refractive index mismatch, 0-matched index
       cfg.isrefint:   1-ref. index mismatch at inner boundaries, [0]-matched index
       cfg.isnormalized:[1]-normalize the output fluence to unitary source, 0-no reflection
       cfg.maxgate:    the num of time-gates per simulation
       cfg.minenergy:  terminate photon when weight less than this level (float) [0.0]
       cfg.unitinmm:   defines the length unit for a grid edge length [1.0]
       cfg.shapes:     a JSON string for additional shapes in the grid
       cfg.internalsrc:set to 1 to skip entry search to speedup launch, effective on AMD GPUs
 
 == GPU settings ==
       cfg.autopilot:  1-automatically set threads and blocks, [0]-use nthread/nblocksize
       cfg.nblocksize: how many CUDA thread blocks to be used [64]
       cfg.nthread:    the total CUDA thread number [2048]
       cfg.gpuid:      which GPU to use (run 'mcx -L' to list all GPUs) [1]
                       if set to an integer, gpuid specifies the index (starts at 1)
                       of the GPU for the simulation; if set to a binary string made
                       of 1s and 0s, it enables multiple GPUs. For example, '1101'
                       allows to use the 1st, 2nd and 4th GPUs together.
       cfg.workload    an array denoting the relative loads of each selected GPU. 
                       for example, [50,20,30] allocates 50%, 20% and 30% photons to the
                       3 selected GPUs, respectively; [10,10] evenly divides the load 
                       between 2 active GPUs. A simple load balancing strategy is to 
                       use the GPU core counts as the weight.
       cfg.isgpuinfo:  1-print GPU info, [0]-do not print
       cfg.sradius:    radius within which we use atomic operations (in grid) [0.0]
                       sradius=0 to disable atomic operations; if sradius=-1,
                       use cfg.crop0 and crop1 to define a cubic atomic zone; if
                       sradius=-2, perform atomic operations in the entire domain;
                       by default, srandius=-2 (atomic operations is used).
 
 == Source-detector parameters ==
       cfg.detpos:     an N by 4 array, each row specifying a detector: [x,y,z,radius]
       cfg.maxdetphoton:   maximum number of photons saved by the detectors [1000000]
       cfg.srctype:    source type, the parameters of the src are specified by cfg.srcparam{1,2}
                       'pencil' - default, pencil beam, no param needed
                       'isotropic' - isotropic source, no param needed
                       'cone' - uniform cone beam, srcparam1(1) is the half-angle in radian
                       'gaussian' [*] - a collimated gaussian beam, srcparam1(1) specifies the waist radius (in voxels)
                       'planar' [*] - a 3D quadrilateral uniform planar source, with three corners specified 
                                 by srcpos, srcpos+srcparam1(1:3) and srcpos+srcparam2(1:3)
                       'pattern' [*] - a 3D quadrilateral pattern illumination, same as above, except
                                 srcparam1(4) and srcparam2(4) specify the pattern array x/y dimensions,
                                 and srcpattern is a floating-point pattern array, with values between [0-1]. 
                       'pattern3d' [*] - a 3D illumination pattern. srcparam1{x,y,z} defines the dimensions,
                                 and srcpattern is a floating-point pattern array, with values between [0-1]. 
                       'fourier' [*] - spatial frequency domain source, similar to 'planar', except
                                 the integer parts of srcparam1(4) and srcparam2(4) represent
                                 the x/y frequencies; the fraction part of srcparam1(4) multiplies
                                 2*pi represents the phase shift (phi0); 1.0 minus the fraction part of
                                 srcparam2(4) is the modulation depth (M). Put in equations:
                                     S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
                       'arcsine' - similar to isotropic, except the zenith angle is uniform
                                 distribution, rather than a sine distribution.
                       'disk' [*] - a uniform disk source pointing along srcdir; the radius is 
                                set by srcparam1(1) (in grid unit)
                       'fourierx' [*] - a general Fourier source, the parameters are 
                                srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phi0,M]
                                normalized vectors satisfy: srcdir cross v1=v2
                                the phase shift is phi0*2*pi
                       'fourierx2d' - a general 2D Fourier basis, parameters
                                srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phix,phiy]
                                the phase shift is phi{x,y}*2*pi
                       'zgaussian' - an angular gaussian beam, srcparam1(0) specifies the variance in the zenith angle
                       'line' - a line source, emitting from the line segment between 
                                cfg.srcpos and cfg.srcpos+cfg.srcparam(1:3), radiating 
                                uniformly in the perpendicular direction
                       'slit' [*] - a colimated slit beam emitting from the line segment between 
                                cfg.srcpos and cfg.srcpos+cfg.srcparam(1:3), with the initial  
                                dir specified by cfg.srcdir
                       'pencilarray' - a rectangular array of pencil beams. The srcparam1 and srcparam2
                                are defined similarly to 'fourier', except that srcparam1(4) and srcparam2(4)
                                are both integers, denoting the element counts in the x/y dimensions, respectively. 
                                For exp., srcparam1=[10 0 0 4] and srcparam2[0 20 0 5] represent a 4x5 pencil beam array
                                spanning 10 grids in the x-axis and 20 grids in the y-axis (5-voxel spacing)
                       source types marked with [*] can be focused using the
                       focal length parameter (4th element of cfg.srcdir)
       cfg.{srcparam1,srcparam2}: 1x4 vectors, see cfg.srctype for details
       cfg.srcpattern: see cfg.srctype for details
       cfg.issrcfrom0: 1-first voxel is [0 0 0], [0]- first voxel is [1 1 1]
       cfg.voidtime:   for wide-field sources, [1]-start timer at launch, or 0-when entering 
                       the first non-zero voxel
 
 == Output control ==
       cfg.issaveexit: [0]-save the position (x,y,z) and (vx,vy,vz) for a detected photon
                       Example: <demo_lambertian_exit_angle.m>
       cfg.issaveref:  [0]-save diffuse reflectance/transmittance in the non-zero voxels
                       next to a boundary voxel. The reflectance data are stored as 
                       negative values; must pad zeros next to boundaries
                       Example: see the demo script at the bottom
       cfg.outputtype: 'flux' - fluence-rate, (default value)
                       'fluence' - fluence integrated over each time gate, 
                       'energy' - energy deposit per voxel
                       'jacobian' or 'wl' - mua Jacobian (replay mode), 
                       'nscat' or 'wp' - weighted scattering counts for computing Jacobian for mus (replay mode)
                       for type jacobian/wl/wp, example: <demo_mcxlab_replay.m>
                       and  <demo_replay_timedomain.m>
       cfg.session:    a string for output file names (only used when no return variables)
 
 == Debug ==
       cfg.debuglevel:  debug flag string, one or a combination of ['R','M','P'], no space
                     'R':  debug RNG, output fluence.data is filled with 0-1 random numbers
                     'M':  return photon trajectory data as the 5th output
                     'P':  show progress bar
       cfg.maxjumpdebug: [1000000|int] when trajectory is requested in the output, 
                      use this parameter to set the maximum position stored. By default,
                      only the first 1e6 positions are stored.
 
       fields with * are required; options in [] are the default values
 
  Output:
       fluence: a struct array, with a length equals to that of cfg.
             For each element of fluence, fluence(i).data is a 4D array with
             dimensions specified by [size(vol) total-time-gates]. 
             The content of the array is the normalized fluence at 
             each voxel of each time-gate.
       detphoton: (optional) a struct array, with a length equals to that of cfg.
             Starting from v2018, the detphoton contains the below subfields:
               detphoton.detid: the ID(>0) of the detector that captures the photon
               detphoton.nscat: cummulative scattering event counts
               detphoton.ppath: cummulative path lengths in each medium (partial pathlength)
                    one need to multiply cfg.unitinmm with ppath to convert it to mm.
               detphoton.p or .v: exit position and direction, when cfg.issaveexit=1
               detphoton.data: a concatenated and transposed array in the order of
                     [detid nscat ppath p v]'
               "data" is the is the only subfield in all MCXLAB-CL before 2018
       vol: (optional) a struct array, each element is a preprocessed volume
             corresponding to each instance of cfg. Each volume is a 3D int32 array.
 
  Example:
       % first query if you have supported GPU(s)
       info=mcxlabcl('gpuinfo')
 
       % define the simulation using a struct
       cfg.nphoton=1e7;
       cfg.vol=uint8(ones(60,60,60));
       cfg.vol(20:40,20:40,10:30)=2;    % add an inclusion
       cfg.prop=[0 0 1 1;0.005 1 0 1.37; 0.2 10 0.9 1.37]; % [mua,mus,g,n]
       cfg.issrcfrom0=1;
       cfg.srcpos=[30 30 1];
       cfg.srcdir=[0 0 1];
       cfg.detpos=[30 20 1 1;30 40 1 1;20 30 1 1;40 30 1 1];
       cfg.vol(:,:,1)=0;   % pad a layer of 0s to get diffuse reflectance
       cfg.issaveref=1;
       cfg.gpuid=1;
       cfg.autopilot=1;
       cfg.tstart=0;
       cfg.tend=5e-9;
       cfg.tstep=5e-10;
       % calculate the fluence distribution with the given config
       [fluence,detpt,vol]=mcxlabcl(cfg);
 
       %%alternatively, you can call
       %[fluence,detpt,vol]=mcxlab(cfg,'mcxcl');
 
       % integrate time-axis (4th dimension) to get CW solutions
       cwfluence=sum(fluence.data,4);  % convert fluence rate to fluence
       cwdref=sum(fluence.dref,4);     % diffuse reflectance
       % plot configuration and results
       subplot(231);
       mcxpreview(cfg);title('domain preview');
       subplot(232);
       mcxplotvol(log10(cwfluence));title('fluence at y=30');
       subplot(233);
       hist(detpt.ppath(:,1),50); title('partial path tissue#1');
       subplot(234);
       plot(squeeze(fluence.data(30,30,30,:)),'-o');title('TPSF at [30,30,30]');
       subplot(235);
       imagesc(squeeze(log(cwdref(:,:,1))));title('diffuse refle. at z=1');
 
  This function is part of Monte Carlo eXtreme (MCX) URL: http://mcx.space/mcxcl
 
  License: GNU General Public License version 3, please read LICENSE.txt for details
</pre>

== # Examples ==

We provided several examples to demonstrate the basic usage of MCXLAB,
as well as to perform validations of MCX algorithm using both simple 
homogeneous and heterogeneous domains. These examples are explained below:

==== demo_mcxlab_basic.m ====

In this example, we show the most basic usage of MCXLAB. This include
how to define the input configuration structure, launch MCX simulations
and interpret and plotting the resulting data.

==== demo_validation_homogeneous.m ====

In this example, we validate MCXLAB with a homogeneous medium in a 
cubic domain. This is exactly the example shown in Fig.5 of [Fang2009].

You can also use the alternative optical properties that has a high g
value to observe the similarity between the two scattering/g configurations.

==== demo_validation_heterogeneous.m ====

In this example, we validate the MCXLAB solver with a heterogeneous
domain and the analytical solution of the diffusion model. We also 
demonstrate how to use sub-pixel resolution to refine the representation
of heterogeneities. The domain is consisted of a 6x6x6 cm box with a 
2cm diameter sphere embedded at the center. 

This test is identical to that used for Fig. 3 in [Fang2010].

==== demo_fullhead_atlas.m ====
In this example, we demonstrate light transport simulation in a full-head 
atlas template(USC 19.5 year group[Sanchez2012]). 
This demo is identical to the MCX simulation used for Fig.9(a) in
[TranYan2019](submitted).

==== demo_mcxyz_skinvessel.m ====
In this example, we compare between MCX and mcxyz written by Dr. Steve Jacques.
The same benchmark can be found at https://omlc.org/software/mc/mcxyz/index.html

==== demo_digimouse_sfdi.m ====
This simulates a widefield SFDI source using the Digimouse atlas. There are
21 tissue types in the atlas.

==== demo_4layer_head.m ====

In this example, we simulate a 4-layer brain model using MCXLAB.
We will investigate the differences between the solutions with and 
witout boundary reflections (both external and internal) and show
you how to display and analyze the resulting data.

==== demo_mcxlab_srctype.m ====
This demo script shows how to use 9 different types sources in your
simulations. These 9 source types include pencil beam, isotropic source,
Gaussian beam, uniform plannar source, uniform disk source, Fourier 
pattern illumuniation (spatial frequency domain sources), arcsine 
distribution beam, uniform cone beam, and an arbitrary light pattern 
(defined by a 2D image).

==== demo_sphere_cube_subpixel.m ====

In this example, we demonstrate how to use sub-pixel resolution 
to represent the problem domain. The domain is consisted of a 
6x6x6 cm box with a 2cm diameter sphere embedded at the center.



== # How to compile MCXLAB-CL ==

To compile MCXLAB-CL for Matlab, you need to cd mcx/src directory, and type 

 make mex

from a shell window. You need to make sure your Matlab is installed and 
the command <tt>mex</tt> is included in your PATH environment variable. Similarly, 
to compile MCXLAB-CL for Octave, you type

 make oct

The command <tt>mkoctfile</tt> must be accessible from your command line
and it is provided in a package named "octave3.x-headers" in Ubuntu (3.x
can be 3.2 or 3.4 etc).

If your graphics card is a Fermi-class or newer, you can compile MCXLAB-CL
with make fermimex or fermioct. The output mex file can determine the
level of atomic operations using the cfg.sradius settings.

== # Screenshots ==

Screenshot for using MCXLAB-CL in Matlab:
  http://mcx.sourceforge.net/upload/matlab_mcxlab.png

Screenshot for using MCXLAB-CL in GNU Octave:
  http://mcx.sourceforge.net/upload/octave_mcxlab.png


== # Reference ==

 [Yu2018] Leiming Yu, Fanny Nina-Paravecino, David Kaeli, Qianqian Fang, 
 "Scalable and massively parallel Monte Carlo photon transport simulations 
 for heterogeneous computing platforms," J. Biomed. Opt. 23(1), 010504 (2018).

 [Fang2009] Qianqian Fang and David A. Boas, "Monte Carlo simulation 
  of photon migration in 3D turbid media accelerated by graphics processing 
  units," Opt. Express 17, 20178-20190 (2009)



