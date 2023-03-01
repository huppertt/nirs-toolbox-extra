% HLMFULLBORNSB Solve full-Born forward problem, infinite slab boundaries
%
%   Calculates the perturbation matrix for the Nth Born Approximation
%
%   [Phi0, A] = HlmFullBornSB(SD, Medium, MeasList, OptProp, Debug)
%
%   A           The forward matrix relating the contribution from each
%               detector or voxel to a specific voxel or detector.  
%                  size(A) = [ (nVox+nDet), (M*nVox+nDet) ];
%
%   Phi0        The incident response at each voxel or detector
%               from each source or detector.  
%                  size(Phi0) = [ (nVox+nDet), (nSrc+nDet) ];
%
%   SD, Medium  The PMI Model structures
%   MeasList    The measurement list for this simulation. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang,
%                     Jonathan Stott
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Phi0, A] = HlmFullBornSB(SD, Medium, MeasList, OptProp, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

% ----------------------------------------------------------------------

%%
%% Calculate local velocity
%%

V = 2.99702458e10 ./ Medium.idxRefr;

%%
%%  Get the extrapolated boundary distances
%%

zBnd = calcExtBnd(Medium.idxRefr, Medium.Muspo);

%%
%% Compute complex wavevectors and diffusion coefficients
%%

if (length(V) ~= length(SD.Lambda))
   error('V and SD.Lambda do not have the same size');
end

D = V ./ (3 * (Medium.Muspo + Medium.Muao));

for f = 1:length(SD.ModFreq)
   % Forward
   K(f,:) = sqrt(V .* Medium.Muao ./ D ...
		 - i * (2*pi * (SD.ModFreq(f) * 1e6)) ./ D);
   % Adjoint
   Kr(f,:)= sqrt(V .* Medium.Muao ./ D ...
		 + i * (2*pi * (SD.ModFreq(f) * 1e6)) ./ D);
end

%%
%% Build the full-born tables
%%

FB = MLtoFB(SD, Medium, MeasList, OptProp);

%%
%% Extract the slab thickness
%%

if (isfield(Medium,'Slab_Thickness'))
   Thickness = Medium.Slab_Thickness;
else
   Thickness = 1;
end

% ----------------------------------------------------------------------

%%
%% Calculate the incident field as the field from each unique source
%% to each unique voxel.  Voxels are implicitly functions of wavelength
%% and frequency because scattering can't mix sources with different 
%% properties.
%%

if (Debug)
   disp('# Incident field');
end

%%
%% Compute the incident response at the voxel/detector from each source.
%%

if (Debug)
   disp('Allocating storage for Phi0');
end

Phi0 = zeros(FB.nPts+FB.nDets, FB.nSrcs+FB.nDets);

for iIdx = 1:FB.nSrcs
   iSrc = FB.MLsrc(iIdx,1);
   iFrq = FB.MLsrc(iIdx,3);
   iWvl = FB.MLsrc(iIdx,4); 
   
   rSrc = moveSrcSlab(SD.SrcPos(iSrc,:), Thickness, Medium.Muspo(iWvl));
   
   if (Debug)
      disp([ 'Source index ' num2str(iIdx) ]);
   end

   % Which set of voxel optical properties does this correspond to?
   iOpt = FB.optR(FB.srcF(iIdx));

   if (FB.MLopt(iOpt,3)~=iFrq | FB.MLopt(iOpt,4)~=iWvl)
      error('MLopt and MLsrc tables out of sync');
   end
   
   % From the optical properties, decide where to put voxel phi
   vvox = (iOpt - 1)*FB.nVox + [1:FB.nVox];

   %% slabA was intended to compute half of the 3-pt function from
   %% optode to voxel, but that's just the 2-pt function really, so I may
   %% as well reused the function.
   
   phi = slabA(rSrc, FB.rVox, 0, D(iWvl), K(iFrq,iWvl), ...
	       V(iWvl), Thickness, zBnd(iWvl)).';

   % Will scale by SD.SrcAmp below
   Phi0(vvox, iIdx) = phi;
   
   % Which set of detectors with these optical properties does this source
   % correspond to?  There's probably a better way to do this using the
   % indices returned by unique(), but this is more obviously correct, so
   % I'll leave it as it is for now.
   
   iMLD = find(FB.MLdet(:,3) == FB.MLopt(iOpt,3) & ...
	       FB.MLdet(:,4) == FB.MLopt(iOpt,4));

   iDet = FB.MLdet(iMLD, 2);
   rDet = moveSrcSlab(SD.DetPos(iDet,:), Thickness, Medium.Muspo(iWvl));
   
   %% If I call the detectors "voxels", slabA will also give me the 2-pt
   %% function from source to detectors, just like I need.
   %%
   %% Rescale these elements by SD.DetAmp since they represent end-points
   
   phi = slabA(rSrc, rDet, 0, D(iWvl), K(iFrq,iWvl), ...
	       V(iWvl), Thickness, zBnd(iWvl)).' .* SD.DetAmp(iDet,iWvl);

   Phi0(FB.nPts+iMLD, iIdx) = phi;
   
   %% Rescale everything by SD.SrcAmp since this is the source term
   Phi0(:,iIdx) = Phi0(:,iIdx) * SD.SrcAmp(iSrc,iWvl);
end

%%
%% Compute the incident response at the voxel from each detector, 
%% these portions are used by Jacobians. 
%%

for iIdx = 1:FB.nDets;
   iDet = FB.MLdet(iIdx,2);
   iFrq = FB.MLdet(iIdx,3);
   iWvl = FB.MLdet(iIdx,4); 
   
   rDet = SD.DetPos(iDet,:);
   
   if (Debug)
      disp([ 'Detector ' num2str(iIdx) ]);
   end
   
   % Which set of voxel optical properties does this correspond to?
   iOpt = FB.optR(FB.detF(iIdx));

   if (FB.MLopt(iOpt,3) ~= iFrq | FB.MLopt(iOpt,4) ~= iWvl)
      error('MLopt and MLsrc tables out of sync');
   end
   
   % From the optical properties, decide where to put voxel phi
   vvox = (iOpt - 1)*FB.nVox + [1:FB.nVox];

   %% slabA was intended to compute half of the 3-pt function from
   %% optode to voxel, but that's just the 2-pt function really, so I may
   %% as well reused the function.
   
   phi = slabA(rDet, FB.rVox, 0, D(iWvl), Kr(iFrq,iWvl), ...
	       V(iWvl), Thickness, zBnd(iWvl)).';
   
   Phi0(vvox,FB.nSrcs + iIdx) = phi;

   % Detector to detector block is all zeros; ignore it
   Phi0(FB.nPts+1:end, FB.nSrcs + iIdx) = 0;
   
   Phi0(:,FB.nSrcs + iIdx) = Phi0(:,FB.nSrcs + iIdx) * SD.DetAmp(iDet,iWvl);
end

% If more than one wavelength/frequency is in use, Phi0 will be sparse
%  with a filling factor roughly nVox/nPts.

if (Debug)
   disp('Sparsifying Phi0');
end

if (FB.nOpt >= 2)
   Phi0 = sparse(Phi0);
end

% ----------------------------------------------------------------------

%% nMeas and nVoxl are not to be taken too literally.

if     (FB.calc_mua & FB.calc_musp)
   nVoxl = 4*FB.nPts;
elseif (FB.calc_musp)
   nVoxl = 3*FB.nPts;
elseif (FB.calc_mua)
   nVoxl =   FB.nPts;
else
   error('OptProp cannot be all zeros');
end

nMeas = FB.nPts + FB.nDets;

% Increased absorption makes the fluence go down
dmua = -1 * FB.volVoxel * ones(length(D),1);

% Want to compute delta-D / D, not delta-D, correct for that.
% Also, actual expression is grad-phi . grad-G, not the
%  grad-phi . grad-phi I computed, so correct for that too.
dmus = D ./ V * FB.volVoxel;

%%
%%  Loop over measurements
%%

if Debug
   disp('Allocating storage for A');
end

% This has to be a dense matrix or everything slows to a crawl while
% filling it.

A = zeros(nMeas, nVoxl+FB.nDets);

if Debug
   disp('Meas #: ');
end

% First, do the measurements that start and end at a voxel.  nPts is 
% nOpt * nVox, so break up the product into two loops for convenience.

for iOpt = 1:FB.nOpt
   iFrq = FB.MLopt(iOpt,3);
   iWvl = FB.MLopt(iOpt,4);
   
   if Debug
      disp([ 'Optical Properties ' num2str(iOpt) ])
   end
   
   for iVox = 1:FB.nVox
      iMeas = (iOpt - 1)*FB.nVox + iVox;
      
      Pos = FB.rVox(iVox,:);
      Amp = 1;
   
      if Debug
	 disp([ 'Voxel ' num2str(iMeas) ])
      end

      %%
      %%  Compute the Green's function from each voxel to each voxel.
      %%

      if (FB.calc_mua)
	 phia = slabA(Pos, FB.rVox, FB.volVoxel, D(iWvl), K(iFrq,iWvl), ...
		      V(iWvl), Thickness, zBnd(iWvl));
      end
   
      if (FB.calc_musp)
	 phis = slabS(Pos, FB.rVox, D(iWvl), K(iFrq,iWvl), ...
		      V(iWvl), Thickness, zBnd(iWvl));
      end
   
      Alst = (iOpt - 1)*FB.nVox + [1:FB.nVox];
   
      if (FB.calc_mua & FB.calc_musp)
	 A(iMeas,             Alst) = phia      * dmua(iWvl);
	 A(iMeas, 1*FB.nPts + Alst) = phis(1,:) * dmus(iWvl);
	 A(iMeas, 2*FB.nPts + Alst) = phis(2,:) * dmus(iWvl);
	 A(iMeas, 3*FB.nPts + Alst) = phis(3,:) * dmus(iWvl);
      elseif (FB.calc_mua)	    		       		 
	 A(iMeas,             Alst) = phia      * dmua(iWvl);
      elseif (FB.calc_musp)	    	       		 
	 A(iMeas,             Alst) = phis(1,:) * dmus(iWvl);
	 A(iMeas, 1*FB.nPts + Alst) = phis(2,:) * dmus(iWvl);
	 A(iMeas, 2*FB.nPts + Alst) = phis(3,:) * dmus(iWvl);
      end
   end
end

% Next, do the measurements that start at a voxel and end at a detector

for iMeas = 1:FB.nDets;
   if Debug
      disp([ 'Detector ' num2str(iMeas) ]);
   end

   iDet = FB.MLdet(iMeas,2);
   iFrq = FB.MLdet(iMeas,3);
   iWvl = FB.MLdet(iMeas,4);

   iOpt = find(FB.MLopt(:,3)==iFrq & FB.MLopt(:,4)==iWvl);
   
   Pos = moveSrcSlab(SD.DetPos(iDet,:), Thickness, Medium.Muspo(iWvl));
   Amp = SD.DetAmp(iDet,iWvl);
   
   %%
   %%  Compute the Green's function from each voxel to each detector.
   %%

   if (FB.calc_mua)
      phia = slabA(Pos, FB.rVox, 0, D(iWvl), Kr(iFrq,iWvl), ...
		   V(iWvl), Thickness, zBnd(iWvl));
   end
   
   if (FB.calc_musp)
      phis = slabS(Pos, FB.rVox, D(iWvl), Kr(iFrq,iWvl), ...
		   V(iWvl), Thickness, zBnd(iWvl));
   end
   
   Alst = (iOpt - 1)*FB.nVox + [1:FB.nVox];
   
   if (FB.calc_mua & FB.calc_musp)
      A(FB.nPts+iMeas,             Alst) = phia      * dmua(iWvl) * Amp;
      A(FB.nPts+iMeas, 1*FB.nPts + Alst) = phis(1,:) * dmus(iWvl) * Amp;
      A(FB.nPts+iMeas, 2*FB.nPts + Alst) = phis(2,:) * dmus(iWvl) * Amp;
      A(FB.nPts+iMeas, 3*FB.nPts + Alst) = phis(3,:) * dmus(iWvl) * Amp;
   elseif (FB.calc_mua)	    		       		 
      A(FB.nPts+iMeas,             Alst) = phia      * dmua(iWvl) * Amp;
   elseif (FB.calc_musp)	    	       		 
      A(FB.nPts+iMeas,             Alst) = phis(1,:) * dmus(iWvl) * Amp;
      A(FB.nPts+iMeas, 1*FB.nPts + Alst) = phis(2,:) * dmus(iWvl) * Amp;
      A(FB.nPts+iMeas, 2*FB.nPts + Alst) = phis(3,:) * dmus(iWvl) * Amp;
   end
end

% Next, do measurements that start at a detector and "end" at a voxel.
% These can be copied out of Phi0.  The have no real physical meaning in
% the forward problem, but I'm not sure about Jacobians.  In any case,
% since they're already computed, the effort is minimal.

A(1:FB.nOpt*FB.nVox, nVoxl + [1:FB.nDets]) = ...
    Phi0(1:FB.nOpt*FB.nVox, FB.nSrcs + [1:FB.nDets]);

% Detector to detector have no physical meaning, so no need to compute
% those terms.

if Debug
   fprintf('\n');
end

% If more than one wavelength/frequency is in use, A will be sparse

if (FB.nOpt >= 2)
   if (Debug)
      disp('Making A sparse');
   end
   
   A = sparse(A);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the incident response from source to detector
%%

function[Phi] = slabA(Src, Vox, volVoxel, D, K, v, Thickness, zBnd)

dR = sum((Vox(:,1)-Src(1)).^2 + ...
	 (Vox(:,2)-Src(2)).^2 + (Vox(:,3)-Src(3)).^2,2);

Phi = zeros(1,length(dR));

nzel = find(dR  > 0);
Vox  = Vox(nzel,:);

% Original source charge.

PhiSrc = FD_GF(Src, Vox, D, K, v); 

% Use a dynamic number of image charges.  Pre-compute as many things
% as possible outside the image-charge loop.
	 
if (Thickness > 0)
   Z0 = 0         - zBnd;
   Z1 = Thickness + zBnd;
else
   Z0 = Thickness - zBnd;
   Z1 = 0         + zBnd;
end

Img = Src;
zs1 = Src(3);
zs2 = Src(3);

bflag = 0;
chg   = 1;

PhiImg = zeros(size(PhiSrc));

for loop = 1:50;               % max 50 image charges
   zi1 = getImageCharge(zs2, Z0);
   zi2 = getImageCharge(zs1, Z1);
	    
   Img(3) = zi1;
   Phi1 = FD_GF(Img, Vox, D, K, v); 

   Img(3) = zi2;
   Phi2 = FD_GF(Img, Vox, D, K, v); 

   if (max(abs(Phi1 + Phi2) ./ abs(PhiSrc + PhiImg)) < 1e-5)
      % Converged sufficiently
      bflag = 1;
   end

   if (abs(Phi1) == 0 & abs(Phi2) == 0)
      % Remaining charges are identically zero
      bflag = 1;
   end
      
   chg  = -chg;

   PhiImg = PhiImg + chg * (Phi1 + Phi2);
	    
   if (bflag ~= 0)
      break;
   else
      % Set up for next iteration
	       
      zs1 = zi1;
      zs2 = zi2;
   end
end

if (bflag == 0)
   error('HlmFullBorn_slab [A] failed to converge');
end

% FD_GF returns a photon density, sum charges and turn into a fluence
Phi(nzel) = v * (PhiSrc + PhiImg).';

% Calculate the self-voxel term (turn voxel into sphere and solve
% analytically).  I'm not conviced the constant factors other constant
% factors are right either.

zel = find(dR == 0);

if length(volVoxel)==1
   R = (volVoxel / (4*pi/3))^(1/3);
else
   R = (volVoxel(zel) / (4*pi/3))^(1/3);
end

% (v/D) or (v/D)^2 ???
Phi(zel) = v/D * (1 - (1 + K*R)*exp(-K*R)) / K^2;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the gradient of response from source to detector
%%

function[GPhi] = slabS(Src, Vox, D, K, v, Thickness, zBnd)

dR = sum((Vox(:,1)-Src(1)).^2 + ...
	 (Vox(:,2)-Src(2)).^2 + (Vox(:,3)-Src(3)).^2,2);

GPhi = zeros(3,length(dR));

nzel = find(dR  > 0);
Vox  = Vox(nzel,:);

% Original source charge (well, gradient there-of).

PhiSrc = FD_GradGF(Src, Vox, D, K, v); 

% Use a dynamic number of image charges.  Pre-compute as many things
% as possible outside the image-charge loop.
	 
if (Thickness > 0)
   Z0 = 0         - zBnd;
   Z1 = Thickness + zBnd;
else
   Z0 = Thickness - zBnd;
   Z1 = 0         + zBnd;
end

Img = Src;
zs1 = Src(3);
zs2 = Src(3);

bflag = 0;
chg   = 1;

PhiImg = zeros(size(PhiSrc));

for loop = 1:50;               % max 50 image charges
   zi1 = getImageCharge(zs2, Z0);
   zi2 = getImageCharge(zs1, Z1);
	    
   Img(3) = zi1;
   Phi1 = FD_GradGF(Img, Vox, D, K, v); 

   Img(3) = zi2;
   Phi2 = FD_GradGF(Img, Vox, D, K, v); 

   dgrad = sqrt( sum(abs(Phi1 + Phi2).^2,2) ...
		 ./ sum(abs(PhiSrc + PhiImg).^2,2));
   
   if (max(dgrad) < 1e-5)
      % Converged sufficiently
      bflag = 1;
   end

   if (all(abs(Phi1) == 0) & all(abs(Phi2) == 0))
      % Remaining charges are identically zero
      bflag = 1;
   end
      
   chg  = -chg;

   PhiImg = PhiImg + chg * (Phi1 + Phi2);
	    
   if (bflag ~= 0)
      break;
   else
      % Set up for next iteration
	       
      zs1 = zi1;
      zs2 = zi2;
   end
end
	 
if (bflag == 0)
   error('Hlm3ptSlab [A] failed to converge');
end

GPhi(:,nzel) = PhiSrc.' + PhiImg.';

% FD_GradGF gave me gradient of photon density, turn it into a fluence
GPhi = v * GPhi;

return;

