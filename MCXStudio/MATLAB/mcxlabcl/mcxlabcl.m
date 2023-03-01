function varargout=mcxlabcl(varargin)
%
%====================================================================
%      MCXLAB-CL - Monte Carlo eXtreme (MCX) for MATLAB/GNU Octave
%--------------------------------------------------------------------
%       Copyright (c) 2018-2019 Qianqian Fang <q.fang at neu.edu>
%                      URL: http://mcx.space
%====================================================================
%
% Format:
%    fluence=mcxlabcl(cfg);
%       or
%    [fluence,detphoton,vol]=mcxlabcl(cfg);
%    [fluence,detphoton,vol]=mcxlabcl(cfg, option);
%
% Input:
%    cfg: a struct, or struct array. Each element of cfg defines 
%         the parameters associated with a simulation. 
%         if cfg='gpuinfo': return the supported GPUs and their parameters,
%         see sample script at the bottom
%    option: (optional), options is a string, specifying additional options
%         option='preview': this plots the domain configuration using mcxpreview(cfg)
%         option='opencl':  force using mcxcl.mex* instead of mcx.mex* on NVIDIA/AMD/Intel hardware
%         option='cuda':    force using mcx.mex* instead of mcxcl.mex* on NVIDIA GPUs
%
%    if one defines USE_MCXCL=1 in MATLAB command line window, all following
%    mcxlab and mcxlabcl calls will use mcxcl.mex; by setting option='cuda', one can
%    force both mcxlab and mcxlabcl to use mcx (cuda version). Similarly, if
%    USE_MCXCL=0, all mcxlabcl and mcxlab call will use mcx.mex by default, unless
%    one set option='opencl'. Please run eval('base','clear USE_MCXCL') to reset.
%
%    cfg may contain the following fields:
%
%== Required ==
%     *cfg.nphoton:    the total number of photons to be simulated (integer)
%     *cfg.vol:        a 3D array specifying the media index in the domain
%                      can be uint8, uint16, uint32, single or double
%                      arrays.
%     *cfg.prop:       an N by 4 array, each row specifies [mua, mus, g, n] in order.
%                      the first row corresponds to medium type 0
%                      (background) which is typically [0 0 1 1]. The
%                      second row is type 1, and so on. The background
%                      medium (type 0) has special meanings: a photon
%                      terminates when moving from a non-zero to zero voxel.
%     *cfg.tstart:     starting time of the simulation (in seconds)
%     *cfg.tstep:      time-gate width of the simulation (in seconds)
%     *cfg.tend:       ending time of the simulation (in second)
%     *cfg.srcpos:     a 1 by 3 vector, the position of the source in grid unit
%     *cfg.srcdir:     a 1 by 3 vector, specifying the incident vector; if srcdir
%                      contains a 4th element, it specifies the focal length of
%                      the source (only valid for focuable src, such as planar, disk,
%                      fourier, gaussian, pattern, slit, etc); if the focal length
%                      is nan, all photons will be launched isotropically regardless
%
%== MC simulation settings ==
%      cfg.seed:       seed for the random number generator (integer) [0]
%                      if set to a uint8 array, the binary data in each column is used 
%                      to seed a photon (i.e. the "replay" mode)
%      cfg.respin:     repeat simulation for the given time (integer) [1]
%      cfg.isreflect:  [1]-consider refractive index mismatch, 0-matched index
%      cfg.isrefint:   1-ref. index mismatch at inner boundaries, [0]-matched index
%      cfg.isnormalized:[1]-normalize the output fluence to unitary source, 0-no reflection
%      cfg.maxgate:    the num of time-gates per simulation
%      cfg.minenergy:  terminate photon when weight less than this level (float) [0.0]
%      cfg.unitinmm:   defines the length unit for a grid edge length [1.0]
%      cfg.shapes:     a JSON string for additional shapes in the grid
%      cfg.internalsrc:set to 1 to skip entry search to speedup launch, effective on AMD GPUs
%
%== GPU settings ==
%      cfg.autopilot:  1-automatically set threads and blocks, [0]-use nthread/nblocksize
%      cfg.nblocksize: how many CUDA thread blocks to be used [64]
%      cfg.nthread:    the total CUDA thread number [2048]
%      cfg.gpuid:      which GPU to use (run 'mcx -L' to list all GPUs) [1]
%                      if set to an integer, gpuid specifies the index (starts at 1)
%                      of the GPU for the simulation; if set to a binary string made
%                      of 1s and 0s, it enables multiple GPUs. For example, '1101'
%                      allows to use the 1st, 2nd and 4th GPUs together.
%      cfg.workload    an array denoting the relative loads of each selected GPU. 
%                      for example, [50,20,30] allocates 50%, 20% and 30% photons to the
%                      3 selected GPUs, respectively; [10,10] evenly divides the load 
%                      between 2 active GPUs. A simple load balancing strategy is to 
%                      use the GPU core counts as the weight.
%      cfg.isgpuinfo:  1-print GPU info, [0]-do not print
%      cfg.optlevel:   optimization level; each level contains all the optimization in the
%                      lower levels; optlevel=2,3 is driver dependent, and may cause dramatic slow down
%                      please avoid or use with caution. Default value is 1.
%                      0 - no optimization; 1- use native math functions; 2/3 - more aggressive optimizations (driver dependent)
%
%== Source-detector parameters ==
%      cfg.detpos:     an N by 4 array, each row specifying a detector: [x,y,z,radius]
%      cfg.maxdetphoton:   maximum number of photons saved by the detectors [1000000]
%      cfg.srctype:    source type, the parameters of the src are specified by cfg.srcparam{1,2}
%                      'pencil' - default, pencil beam, no param needed
%                      'isotropic' - isotropic source, no param needed
%                      'cone' - uniform cone beam, srcparam1(1) is the half-angle in radian
%                      'gaussian' [*] - a collimated gaussian beam, srcparam1(1) specifies the waist radius (in voxels)
%                      'planar' [*] - a 3D quadrilateral uniform planar source, with three corners specified 
%                                by srcpos, srcpos+srcparam1(1:3) and srcpos+srcparam2(1:3)
%                      'pattern' [*] - a 3D quadrilateral pattern illumination, same as above, except
%                                srcparam1(4) and srcparam2(4) specify the pattern array x/y dimensions,
%                                and srcpattern is a floating-point pattern array, with values between [0-1]. 
%                      'pattern3d' [*] - a 3D illumination pattern. srcparam1{x,y,z} defines the dimensions,
%                                and srcpattern is a floating-point pattern array, with values between [0-1]. 
%                      'fourier' [*] - spatial frequency domain source, similar to 'planar', except
%                                the integer parts of srcparam1(4) and srcparam2(4) represent
%                                the x/y frequencies; the fraction part of srcparam1(4) multiplies
%                                2*pi represents the phase shift (phi0); 1.0 minus the fraction part of
%                                srcparam2(4) is the modulation depth (M). Put in equations:
%                                    S=0.5*[1+M*cos(2*pi*(fx*x+fy*y)+phi0)], (0<=x,y,M<=1)
%                      'arcsine' - similar to isotropic, except the zenith angle is uniform
%                                distribution, rather than a sine distribution.
%                      'disk' [*] - a uniform disk source pointing along srcdir; the radius is 
%                               set by srcparam1(1) (in grid unit)
%                      'fourierx' [*] - a general Fourier source, the parameters are 
%                               srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phi0,M]
%                               normalized vectors satisfy: srcdir cross v1=v2
%                               the phase shift is phi0*2*pi
%                      'fourierx2d' - a general 2D Fourier basis, parameters
%                               srcparam1: [v1x,v1y,v1z,|v2|], srcparam2: [kx,ky,phix,phiy]
%                               the phase shift is phi{x,y}*2*pi
%                      'zgaussian' - an angular gaussian beam, srcparam1(0) specifies the variance in the zenith angle
%                      'line' - a line source, emitting from the line segment between 
%                               cfg.srcpos and cfg.srcpos+cfg.srcparam(1:3), radiating 
%                               uniformly in the perpendicular direction
%                      'slit' [*] - a colimated slit beam emitting from the line segment between 
%                               cfg.srcpos and cfg.srcpos+cfg.srcparam(1:3), with the initial  
%                               dir specified by cfg.srcdir
%                      'pencilarray' - a rectangular array of pencil beams. The srcparam1 and srcparam2
%                               are defined similarly to 'fourier', except that srcparam1(4) and srcparam2(4)
%                               are both integers, denoting the element counts in the x/y dimensions, respectively. 
%                               For exp., srcparam1=[10 0 0 4] and srcparam2[0 20 0 5] represent a 4x5 pencil beam array
%                               spanning 10 grids in the x-axis and 20 grids in the y-axis (5-voxel spacing)
%                      source types marked with [*] can be focused using the
%                      focal length parameter (4th element of cfg.srcdir)
%      cfg.{srcparam1,srcparam2}: 1x4 vectors, see cfg.srctype for details
%      cfg.srcpattern: see cfg.srctype for details
%      cfg.issrcfrom0: 1-first voxel is [0 0 0], [0]- first voxel is [1 1 1]
%      cfg.voidtime:   for wide-field sources, [1]-start timer at launch, or 0-when entering 
%                      the first non-zero voxel
%
%== Output control ==
%      cfg.issaveexit: [0]-save the position (x,y,z) and (vx,vy,vz) for a detected photon
%                      Example: <demo_lambertian_exit_angle.m>
%      cfg.issaveref:  [0]-save diffuse reflectance/transmittance in the non-zero voxels
%                      next to a boundary voxel. The reflectance data are stored as 
%                      negative values; must pad zeros next to boundaries
%                      Example: see the demo script at the bottom
%      cfg.outputtype: 'flux' - fluence-rate, (default value)
%                      'fluence' - fluence integrated over each time gate, 
%                      'energy' - energy deposit per voxel
%                      'jacobian' or 'wl' - mua Jacobian (replay mode), 
%                      'nscat' or 'wp' - weighted scattering counts for computing Jacobian for mus (replay mode)
%                      for type jacobian/wl/wp, example: <demo_mcxlab_replay.m>
%                      and  <demo_replay_timedomain.m>
%      cfg.session:    a string for output file names (only used when no return variables)
%
%== Debug ==
%      cfg.debuglevel:  debug flag string, one or a combination of ['R','M','P'], no space
%                    'R':  debug RNG, output fluence.data is filled with 0-1 random numbers
%                    'M':  return photon trajectory data as the 5th output
%                    'P':  show progress bar
%      cfg.maxjumpdebug: [1000000|int] when trajectory is requested in the output, 
%                     use this parameter to set the maximum position stored. By default,
%                     only the first 1e6 positions are stored.
%
%      fields with * are required; options in [] are the default values
%
% Output:
%      fluence: a struct array, with a length equals to that of cfg.
%            For each element of fluence, fluence(i).data is a 4D array with
%            dimensions specified by [size(vol) total-time-gates]. 
%            The content of the array is the normalized fluence at 
%            each voxel of each time-gate.
%      detphoton: (optional) a struct array, with a length equals to that of cfg.
%            Starting from v2018, the detphoton contains the below subfields:
%              detphoton.detid: the ID(>0) of the detector that captures the photon
%              detphoton.nscat: cummulative scattering event counts
%              detphoton.ppath: cummulative path lengths in each medium (partial pathlength)
%                   one need to multiply cfg.unitinmm with ppath to convert it to mm.
%              detphoton.p or .v: exit position and direction, when cfg.issaveexit=1
%              detphoton.data: a concatenated and transposed array in the order of
%                    [detid nscat ppath p v]'
%              "data" is the is the only subfield in all MCXLAB-CL before 2018
%      vol: (optional) a struct array, each element is a preprocessed volume
%            corresponding to each instance of cfg. Each volume is a 3D int32 array.
%
% Example:
%      % first query if you have supported GPU(s)
%      info=mcxlab('gpuinfo')
%
%      % define the simulation using a struct
%      cfg.nphoton=1e7;
%      cfg.vol=uint8(ones(60,60,60));
%      cfg.vol(20:40,20:40,10:30)=2;    % add an inclusion
%      cfg.prop=[0 0 1 1;0.005 1 0 1.37; 0.2 10 0.9 1.37]; % [mua,mus,g,n]
%      cfg.issrcfrom0=1;
%      cfg.srcpos=[30 30 1];
%      cfg.srcdir=[0 0 1];
%      cfg.detpos=[30 20 1 1;30 40 1 1;20 30 1 1;40 30 1 1];
%      cfg.vol(:,:,1)=0;   % pad a layer of 0s to get diffuse reflectance
%      cfg.issaveref=1;
%      cfg.gpuid=1;
%      cfg.autopilot=1;
%      cfg.tstart=0;
%      cfg.tend=5e-9;
%      cfg.tstep=5e-10;
%      % calculate the fluence distribution with the given config
%      [fluence,detpt,vol]=mcxlabcl(cfg);
%
%      %%alternatively, you can call
%      %[fluence,detpt,vol]=mcxlab(cfg,'mcxcl');
%
%      % integrate time-axis (4th dimension) to get CW solutions
%      cwfluence=sum(fluence.data,4);  % convert fluence rate to fluence
%      cwdref=sum(fluence.dref,4);     % diffuse reflectance
%      % plot configuration and results
%      subplot(231);
%      mcxpreview(cfg);title('domain preview');
%      subplot(232);
%      imagesc(squeeze(log(cwfluence(:,30,:))));title('fluence at y=30');
%      subplot(233);
%      hist(detpt.ppath(:,1),50); title('partial path tissue#1');
%      subplot(234);
%      plot(squeeze(fluence.data(30,30,30,:)),'-o');title('TPSF at [30,30,30]');
%      subplot(235);
%      imagesc(squeeze(log(cwdref(:,:,1))));title('diffuse refle. at z=1');
%      mcxplotvol(log10(cwfluence));title('fluence at y=30');
%
% This function is part of Monte Carlo eXtreme (MCX) URL: http://mcx.space/mcxcl
%
% License: GNU General Public License version 3, please read LICENSE.txt for details
%

try
    defaultocl=evalin('base','USE_MCXCL');
catch
    defaultocl=1;
end

useopencl=defaultocl;

if(nargin==2 && ischar(varargin{2}))
    if(strcmp(varargin{2},'preview'))
        [varargout{1:nargout}]=mcxpreview(varargin{1});
	    return;
    elseif(strcmp(varargin{2},'cuda'))
        useopencl=0;
    end
end

if(isstruct(varargin{1}))
    for i=1:length(varargin{1})
        castlist={'srcpattern','srcpos','detpos','prop','workload','srcdir'};
        for j=1:length(castlist)
            if(isfield(varargin{1}(i),castlist{j}))
                varargin{1}(i).(castlist{j})=double(varargin{1}(i).(castlist{j}));
            end
        end
    end
end

if(useopencl==0)
    [varargout{1:nargout}]=mcx(varargin{1});
else
    [varargout{1:nargout}]=mcxcl(varargin{1});
end

if(nargin==0)
    return;
end

cfg=varargin{1};

if(nargout>=2)

    for i=1:length(varargout{2})
        if(~isfield(cfg(i),'issaveexit') || cfg(i).issaveexit~=2)
            medianum=size(cfg(i).prop,1)-1;
            detp=varargout{2}(i).data;
            if(isempty(detp))
                continue;
            end
            newdetp.detid=int32(detp(1,:))';
            newdetp.nscat=int32(detp(2,:))';    % 1st medianum block is num of scattering
            newdetp.ppath=detp(3:2+medianum,:)';% 2nd medianum block is partial path
            if(isfield(cfg(i),'issaveexit') && cfg(i).issaveexit)
                newdetp.p=detp(end-5:end-3,:)';             %columns 7-5 from the right store the exit positions*/
                newdetp.v=detp(end-2:end,:)';	     %columns 4-2 from the right store the exit dirs*/
            end
            % newdetp.w0=detp(end,:)';  % last column is the initial packet weight
            newdetp.data=detp;      % enable this line for compatibility
            newdetpstruct(i)=newdetp;
        else
            newdetpstruct(i)=varargout{2}(i);
        end
    end
    if(exist('newdetpstruct','var'))
        varargout{2}=newdetpstruct;
    end

    if(nargout>=5)
        for i=1:length(varargout{5})
            data=varargout{5}.data;
            if(isempty(data))
               continue;
            end
            traj.pos=data(2:4,:).';
            traj.id=typecast(data(1,:),'uint32').';
            [traj.id,idx]=sort(traj.id);
            traj.pos=traj.pos(idx,:);
            traj.data=[single(traj.id)' ; data(2:end,idx)];
            newtraj(i)=traj;
        end
        if(exist('newtraj','var'))
            varargout{5}=newtraj;
        end
    end
end
