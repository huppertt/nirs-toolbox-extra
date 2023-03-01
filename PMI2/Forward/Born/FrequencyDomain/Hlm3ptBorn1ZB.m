% HLM3PTBORN1ZB  1st Born approx. of a perturbed Helmholtz equation
%                in a semi-infinite geometry
%
% HLM3PTBORN1ZB computes the forward weighting matrix associated
% with the Born-1 approximation to a spatial varying k.  A zero
% fluence boudary condition is implemented at the distance zBnd by
% mirroring the sources and the pertubation responses.  The sources
% are assumed to be unit amplitude point sources.  Absorption
% reconstructs delta-mu, scattering reconstructs delta-D/D.
% 
% [Phi0,A] = Hlm3ptBorn1ZB(SD, Medium, MeasList, OptProp, Debug)
%
%   SD, Medium  The PMI Model structures
%   MeasList    The measurement list for this simulation.
%   OptProp     Perturbations to model ([ Absflag, Scatflag ])
%   Debug       OPTIONAL: Print out debugging info.
%
% Returns:
%   A           The forward matrix relating the contribution from each voxel
%               to a specific source-detector pair.  Each row is for a
%               different source detector pair.
%   Phi0        The incident response at each detector from each source in
%               a column vector.  The same combination pattern as the
%               columns of A.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2004, David Boas, Dana Brooks, Rick Gaudette, 
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

function [Phi0, A] = Hlm3ptBorn1ZB(SD, Medium, MeasList, OptProp, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var'))
   Debug = 0;
end

%%
%% Calculate local velocity
%%

V = 2.99702458e10 ./ Medium.idxRefr;

nMeas = size(MeasList,1);

%%
%%  Get the extrapolated boundary distance
%%

zBnd = calcExtBnd(Medium.idxRefr, Medium.Muspo);

%%
%% Compute complex wavevectors and diffusion coefficients
%%

if (length(V) < max(unique(MeasList(:,4))))
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
%% Extract the slab thickness
%%

Thickness = Medium.Slab_Thickness;

%%
%%  Create the sampling volume
%%

[Xm Ym Zm volVoxel] = sampleVolume(Medium.CompVol);

nVox = length(Xm);
rVox = [ Xm Ym Zm ];    % Nx3 vector of voxel positions
clear Xm Ym Zm;

calc_mua  = OptProp(1);
calc_musp = OptProp(2);

%%
%%  Loop over each source detector combination
%%

vsrc = unique(MeasList(:,1));
vdet = unique(MeasList(:,2));
vfrq = unique(MeasList(:,3));
vwvl = unique(MeasList(:,4));

% --------------------------------------

if (calc_mua)
   if (Debug)
      disp('Calculating absorbing matrix terms');
   end
   
   Amua = zeros(nMeas, nVox);

   % Do all the source 2-pts
   
   for f = 1:length(vfrq)
      iFrq = vfrq(f);
      
      for s = 1:length(vsrc)
	 iSrc = vsrc(s);
	 
	 if (Debug)
	    disp([ 'Source ' num2str(iSrc) ', frequency ' num2str(iFrq) ]);
	 end
	 
	 for l = 1:length(vwvl)
	    iWav = vwvl(l);
	    
	    ml = find(MeasList(:,1)==iSrc & ...
		      MeasList(:,3)==iFrq & MeasList(:,4)==iWav);
	    
	    if (~isempty(ml))
	       % Move the effective optode positions one mean free path 
	       % into the medium.

	       Src = moveSrcSlab(SD.SrcPos(iSrc,:), ...
				 sign(Thickness)*inf, Medium.Muspo(iWav));
	       
	       phiS = slabA(Src, rVox, D(iWav), K(iFrq, iWav), ...
			    V(iWav), Thickness, zBnd(iWav));
	       phiS = reshape(phiS, 1, nVox) * SD.SrcAmp(iSrc, iWav, iFrq);
	       
	       for m = 1:length(ml)
		  Amua(ml(m),:) = phiS;
	       end
	    end
	 end
      end
   end
   
   % Multiply through by all the detector 2-pts
   
   for f = 1:length(vfrq)
      iFrq = vfrq(f);
      
      for k = 1:length(vdet)
	 iDet = vdet(k);
	 
	 if (Debug)
	    disp([ 'Detector ' num2str(iDet) ', frequency ' num2str(iFrq) ]);
	 end
	 
	 for l = 1:length(vwvl)
	    iWav = vwvl(l);
	    
	    ml = find(MeasList(:,2)==iDet & ...
		      MeasList(:,3)==iFrq & MeasList(:,4)==iWav);
	    
	    if (~isempty(ml))
	       % Move the effective optode positions one mean free path 
	       % into the medium.

	       Det = moveSrcSlab(SD.DetPos(iDet,:), ...
				 sign(Thickness)*inf, Medium.Muspo(iWav));
	       
	       phiD = slabA(Det, rVox, D(iWav), Kr(iFrq, iWav), ...
			    V(iWav), Thickness, zBnd(iWav));
	       phiD = reshape(conj(phiD),1,nVox) * SD.DetAmp(iDet, iWav, iFrq);
	       
	       for m = 1:length(ml)
		  Amua(ml(m),:) = Amua(ml(m),:) .* phiD;
	       end
	    end
	 end
      end
   end

   % Increased absorption makes the fluence go down
   Amua = -1 * Amua;

   clear iDet iFrq iSrc iWav f k l m ml phiS phiD s Src Det;
else
   Amua = [];
end

% --------------------------------------

if (calc_musp)
   if (Debug)
      disp('Calculating scattering matrix terms');
   end
   
   Amus = zeros(nMeas, nVox, 3);

   % Do all the source 2-pts
   
   for f = 1:length(vfrq)
      iFrq = vfrq(f);
	    
      for k = 1:length(vsrc)
	 iSrc = vsrc(k);
      
	 if (Debug)
	    disp([ 'Source ' num2str(iSrc) ', frequency ' num2str(iFrq) ]);
	 end
      
	 for l = 1:length(vwvl)
	    iWav = vwvl(l);
	 
	    ml = find(MeasList(:,1)==iSrc & ...
		      MeasList(:,3)==iFrq & MeasList(:,4)==iWav);
	       
	    if (~isempty(ml))
	       % Move the effective optode positions one mean free path 
	       % into the medium.

	       Src = moveSrcSlab(SD.SrcPos(iSrc,:), ...
				 sign(Thickness)*inf, Medium.Muspo(iWav));
	       
	       phiS = slabS(Src, rVox, D(iWav), K(iFrq, iWav), ...
			    V(iWav), Thickness, zBnd(iWav));
	       phiS = reshape(phiS, 1, nVox, 3) * SD.SrcAmp(iSrc, iWav, iFrq);
	       
	       for m = 1:length(ml)
		  Amus(ml(m),:,1:3) = phiS;
	       end
	    end
	 end
      end
   end
   
   % Multiply through by all the detector 2-pts
   
   for f = 1:length(vfrq)
      iFrq = vfrq(f);
      
      for k = 1:length(vdet)
	 iDet = vdet(k);
	 
	 if (Debug)
	    disp([ 'Detector ' num2str(iDet) ', frequency ' num2str(iFrq) ]);
	 end
	 
	 for l = 1:length(vwvl)
	    iWav = vwvl(l);
	    
	    ml = find(MeasList(:,2)==iDet & ...
		      MeasList(:,3)==iFrq & MeasList(:,4)==iWav);
	    
	    if (~isempty(ml))
	       % Move the effective optode positions one mean free path 
	       % into the medium.

	       Det = moveSrcSlab(SD.DetPos(iDet,:), ...
				 Thickness, Medium.Muspo(iWav));
	       
	       phiD = slabS(Det, rVox, D(iWav), Kr(iFrq, iWav), ...
			    V(iWav), Thickness, zBnd(iWav));
	       phiD = reshape(conj(phiD), 1, nVox, 3) * ...
		      SD.DetAmp(iDet, iWav, iFrq);
	       
	       for m = 1:length(ml)
		  Amus(ml(m),:,:) = Amus(ml(m),:,:) .* phiD;
	       end
	    end
	 end
      end
   end

   % Take the dot product

   if (Debug)
      disp('Collapsing vectors');
   end
   
   Amus = sum(Amus,3);
   
   % Actual expression is grad-phi . grad-G, not the
   %  grad-phi . grad-phi I computed, so correct for that here
   
   for m = 1:size(Amus,1)
      iWav = MeasList(m,4);
      
      Amus(m,:) = Amus(m,:) * D(iWav) / V(iWav);
   end

   clear iDet iFrq iSrc iWav k l m ml phiS phiD s f Src Det;
else
   Amus = [];
end

clear vsrc vdet vfrq vwvl;

% --------------------------------------

if (Debug)
   disp('Generating final matrix');
end

A = [ Amua Amus ] * volVoxel;

clear Amua Amus;

% Rather than re-inventing the whell, just call the same 2-pt code as
% everything else.

Phi0 = FD2pt(SD,Medium,MeasList,Debug);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the incident response from source to detector
%%

function[Phi] = slabA(Src, Vox, D, K, v, Thickness, zBnd)

% Original source charge.

PhiSrc = FD_GF(Src, Vox, D, K, v); 

% Image charge

Img = Src;

if (Thickness > 0)
   Img(3) = getImageCharge(Src(3), -zBnd);
else
   Img(3) = getImageCharge(Src(3),  zBnd);
end

PhiImg = -1*FD_GF(Img, Vox, D, K, v); 

% Multiply by charge and sum to get final field

Phi = PhiSrc + PhiImg;

% Turn photon density into fluence
Phi = v * Phi;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the gradient of response from source to detector
%%

function[GPhi] = slabS(Src, Vox, D, K, v, Thickness, zBnd)
	
% Original source charge (well, gradient there-of).  Green's functions
% are for photon density, multiply by v to convert to fluence.

PhiSrc = FD_GradGF(Src, Vox, D, K, v); 

% Image charge

Img = Src;

if (Thickness > 0)
   Img(3) = getImageCharge(Src(3), -zBnd);
else
   Img(3) = getImageCharge(Src(3),  zBnd);
end

PhiImg = FD_GradGF(Img, Vox, D, K, v); 

% Multiply by charge and sum to get final field

GPhi = PhiSrc - PhiImg;

% Turn photon density into fluence
GPhi = v * GPhi;

return;
