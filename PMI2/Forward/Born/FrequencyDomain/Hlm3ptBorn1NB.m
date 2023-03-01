% HLM3PTBORN1NB  1st Born approx. of a perturbed Helmholtz equation
%                with infinite boundaries.
%
% HLM3PTBORN1NB computes the forward weighting matrix associated
% with the Born-1 approximation to a spatial varying k with
% infinte boundary conditions.  The sources are assumed to be unit 
% amplitude point sources.  Absorption reconstructs delta-mu, 
% scattering reconstructs delta-D/D.
% 
% [Phi0,A] = Hlm3ptBorn1NB(SD, Medium, MeasList, OptProp, Debug)
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

function [Phi0, A] = Hlm3ptBorn1NB(SD, Medium, MeasList, OptProp, Debug)

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
%%  Create the sampling volume
%%

[Xm Ym Zm volVoxel] = sampleVolume(Medium.CompVol);

nVox = length(Xm);
rVox = [ Xm Ym Zm ];    % Nx3 vector of voxel positions

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
	       phiS = slabA(SD.SrcPos(iSrc,:), rVox, D(iWav), ...
			    K(iFrq, iWav), V(iWav));
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
	       phiD = slabA(SD.DetPos(iDet,:), rVox, D(iWav), ...
			    Kr(iFrq, iWav), V(iWav));
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

   clear iDet iFrq iSrc iWav k l m ml phiS phiD s f Src Det;
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
	       phiS = slabS(SD.SrcPos(iSrc,:), rVox, D(iWav), ...
			    K(iFrq, iWav), V(iWav));
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
	       phiD = slabS(SD.DetPos(iDet,:), rVox, D(iWav), ...
			    Kr(iFrq, iWav), V(iWav));
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

function[Phi] = slabA(Src, Vox, D, K, v)

% FD_GF returns a photon density, turn it into a fluence

Phi = v * FD_GF(Src, Vox, D, K, v); 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the gradient of response from source to detector
%%

function[GPhi] = slabS(Src, Vox, D, K, v)

% FD_GradGF gave me gradient of photon density, turn it into a fluence
% (well, gradient there-of)

GPhi = v * FD_GradGF(Src, Vox, D, K, v); 

return;


return
