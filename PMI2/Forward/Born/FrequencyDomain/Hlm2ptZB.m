%HLM2PTZB  Exact solution of Helmholtz equation for semi-infinite boundaries.
%
%   Phi0 = Hlm2ptNB(SD, Medium, MeasList)
%
%   Phi0     The incident response, in measurement list order.
%
%   SD, Medim   The PMI Model structures
%   MeasList    The measurement list for this simulation.
%
%   HLM2PTZB computes the exact solution of the Helmholtz equation
%   for a spatial uniform K.  Semi-infinite boundaries, interface
%   implicitly define the Z=0 plane.

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

function[Phi0] = Hlm2ptZB(SD, Medium, MeasList, Debug)

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
   %%  Move sources and detectors one mean free path into medium
   %%

   Src = moveSrcSlab(Src, inf, Medium.Muspo(iWav));
   Det = moveSrcSlab(Det, inf, Medium.Muspo(iWav));
   
   %%
   %%  Compute location of image charge
   %%

   Img    = Src;
   Img(3) = getImageCharge(Src(3), -zBnd(iWav));
   
   %%
   %% Compute the incident response at the detector, rescale using
   %% amplitude factors from SD structure, and convert from photon density
   %% (J/cm^3) to fluence (W/cm^2) by multiplying by the speed of light.
   %%
   
   PhiSrc = FD_GF(Src, Det, D(iWav), K(iFrq, iWav), V(iWav));
   PhiImg = FD_GF(Img, Det, D(iWav), K(iFrq, iWav), V(iWav));

   Phi0(iMeas) = (PhiSrc - PhiImg) * V(iWav) ...
                 * SD.SrcAmp(iSrc,iWav,iFrq) * SD.DetAmp(iDet,iWav,iFrq);
end

if Debug
   fprintf(1, '\n');
end

return;
