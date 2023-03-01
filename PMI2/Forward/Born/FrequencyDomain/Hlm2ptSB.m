% HLM2PTSB  Exact solution of Helmholtz equation for a slab.
%
%   Phi0 = Hlm2ptSB(SD, Medium, MeasList)
%
%   Phi0     The incident response at each detector from each source
%               in a column vector.  The same combination pattern as 
%               the columns of A.
%
%   SD, Medim   The PMI Model structures
%   MeasList    The measurement list for this simulation.
%
%   HLM2PTSLAB computes the exact solution of the Helmholtz equation
%   for a spatial uniform k.  A slab boudary condition is implemented
%   by mirroring the sources.

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

function[Phi0] = Hlm2ptSB(SD, Medium, MeasList, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

%%
%% Calculate local velocity
%%

V = 2.99702458e10 ./ Medium.idxRefr;

nMeas = size(MeasList,1);

%%
%% Extract the slab thickness
%%

if (Medium.Slab_Thickness == 0)
   error('Illegal zero slab thickness detected');
end

Thickness = Medium.Slab_Thickness;

%% Initialize as many variables as possible outside the loop

Phi0 = zeros(nMeas, 1);

D = V ./ (3 * (Medium.Muspo + Medium.Muao));
	 
for ifrq = 1:length(SD.ModFreq)
   K(ifrq,:) = sqrt(V .* Medium.Muao ./ D ...
	       - i * (2*pi * (SD.ModFreq(ifrq) * 1e6)) ./ D);
end

zBnd = calcExtBnd(Medium.idxRefr, Medium.Muspo);

%%
%%  Loop over each source detector combination
%%

if Debug
   fprintf(1, 'Meas #: ');
end

for iMeas = 1:size(MeasList,1)
   iSrc = MeasList(iMeas,1);
   iDet = MeasList(iMeas,2);
   iFrq = MeasList(iMeas,3);
   iWav = MeasList(iMeas,4);

   Src = SD.SrcPos(iSrc,:);
   Det = SD.DetPos(iDet,:);
   
   if Debug
      fprintf(1, '%d  ', iMeas);
   end
   
   %%
   %%  Move the effective source position 1 mean free path into the medium
   %%

   Src = moveSrcSlab(Src, Thickness, Medium.Muspo(iWav));
   Det = moveSrcSlab(Det, Thickness, Medium.Muspo(iWav));
   
   %%
   %%  Compute the incident response at the detector.
   %%
   
   PhiSrc = FD_GF(Src, Det, D(iWav), K(iFrq, iWav), V(iWav)); 

   if (all(PhiSrc(:) == 0))
      % No number of image charges will fix this
      Phi0 = PhiSrc;
      return;
   end
   
   %%
   %% Use a dynamic number of image charges.  Pre-compute as much
   %% as possible outside the image-charge loop.
   %%

   if (Thickness > 0)
      Z0 = 0         - zBnd(iWav);
      Z1 = Thickness + zBnd(iWav);
   else
      Z0 = 0         + zBnd(iWav);
      Z1 = Thickness - zBnd(iWav);
   end
   
   bflag = 0;
   chg   = 1;

   Img = Src;
   
   zs1 = Src(3); % Start at source z-value
   zs2 = Src(3);
   
   PhiImg = 0;
   
   for loop = 1:50;               % max 50 image charges
      zi1 = getImageCharge(zs2, Z0);
      zi2 = getImageCharge(zs1, Z1);
      
      Img(3) = zi1;
      Phi1 = FD_GF(Img, Det, D(iWav), K(iFrq, iWav), V(iWav)); 

      Img(3) = zi2;
      Phi2 = FD_GF(Img, Det, D(iWav), K(iFrq, iWav), V(iWav)); 
      
      if (abs(Phi1 + Phi2) ./ abs(PhiSrc + PhiImg) < 1e-5)
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
      error('Hlm2ptSB failed to converge');
   end

   %%
   %% Compute the incident response at the detector, rescale using
   %% amplitude factors from SD structure, and convert from photon density
   %% (J/cm^3) to fluence (W/cm^2) by multiplying by the speed of light.
   %%
   
   Phi0(iMeas) = (PhiSrc + PhiImg) * V(iWav) ...
                      * SD.SrcAmp(iSrc,iWav,iFrq) * SD.DetAmp(iDet,iWav,iFrq);
end

if Debug
   fprintf(1, '\n');
end

return;
