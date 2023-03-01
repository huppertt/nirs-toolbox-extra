% fitSDPhs   Compute Source-Detector scaling vector, complex phase only
%
% [SrcPhs, DetPhs, DataPhs] = fitSDPhs(SD, Medium, MeasList, data, [Phi0]);
%
% Inputs:
%   SD, Medium - PMI2 structures
%   data       - Measured fluence vector.  The length of the vector
%                must match length(find(SD.MeasList(:,4)==l)).
%   l          - wavelength to fit
%
% Outputs:
%   SrcPhs, DetPhs - per-source and per-detector scaling
%                coefficients.  These can be copied directly into
%                the SD structure, if desired.
%   DataPhs        - per-measurement scaling coefficients, derived
%                from SrcPhs, DetPhs and SD.MeasList.  Multiply
%                the incident fluence computed from the supplied SD
%                by DataPhs to get the re-scaled incident fluence.

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

function[SrcPhs, DetPhs, DataPhs] = fitSDPhs(SD, Medium, MeasList, ...
					     data, Phi0, Debug);

% Copy out of SD if missing or not provided
if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

if (length(unique(MeasList(:,1))) == 1)
   error('SD fitting requires more than one source');
end

if (length(unique(MeasList(:,2))) == 1)
   error('SD fitting requires more than one detector');
end

if (~exist('Phi0','var') | isempty(Phi0))
   if (Debug)
      disp('Calculating forward problem');
   end
   
   % Get Phi0 for current SD values
   Phi0 = DPDWHelmholtz(SD, Medium, MeasList);
end

if (any(size(data) ~= size(Phi0)))
   warning('Lengths of fluence vectors do not match');
   keyboard;
end

SrcPhs  = zeros(size(SD.SrcAmp));
DetPhs  = zeros(size(SD.DetAmp));

% Get Jacobian for computing SD coefficients in the Rytov approx

if (Debug)
   disp('Calculating amplitude Jacobian');
end

Jamp = genSDJacobian(SD, Medium, MeasList, [], 'Rytov');

% The matrix as constructed above, is rank-deficient.  The simplest
% solution is to knock out the first detector used at each wavelength.
% This fails to work, however, when we're combining data sets from two
% imagers (which the program has now way to know about), since the data
% will take on a block form (i.e. the imagers are blind to each other).  I
% don't see any good way around this problem.

ns = size(SD.SrcPos,1);
nd = size(SD.DetPos,1);
nf = max(1,length(SD.ModFreq));
nw = max(1,length(SD.Lambda));

% BUG: If two sub-sets of the optodes do not talk to each other
% (specifically, if the matrix Jamp has a block form) then we will 
% still be rank deficient!

for l = 1:nw
   for m = 1:nf
      il = find(MeasList(:,4)==l & MeasList(:,3)==m);
   
      if ~isempty(il)
	 di = min(MeasList(il,2));
      
	 di = di + (m-1)*nd*nw + (l-1)*nd + nw*ns*nf;
	 Jamp(il, di) = 0;
      end
   end
end

% Eliminate missing (sum==0) and singular (sum==1) rows

nzel = find(sum(Jamp.^2) > 1);

% Compute scattered fluence in Rytov approx, phase only

Phi1 = imag(log(data ./ Phi0));

% Solve the system of linear equations

if (Debug)
   disp('Inverting for optode coefficients');
end

sdfit = zeros(size(Jamp,2),1);
sdfit(nzel) = Jamp(:,nzel) \ Phi1;

% Repack to have the right size

if (Debug)
   disp('Repacking computed coefficients');
end

SrcPhs = reshape(sdfit(         [1:nw*ns*nf]), ns, nw, nf);
DetPhs = reshape(sdfit(nw*nf*ns+[1:nw*nd*nf]), nd, nw, nf);

clear Jamp sdfit

% What was computed is actually the log of the scaling
% coefficients, transform to a linear scale.

SrcPhs = exp(i*SrcPhs);
DetPhs = exp(i*DetPhs);

% Compute the vector for rescaling Phi0

for k = 1:length(MeasList)
   si = MeasList(k,1);
   di = MeasList(k,2);
   fi = max(1,MeasList(k,3));
   wi = max(1,MeasList(k,4));
   
   DataPhs(k) = SrcPhs(si,wi,fi) * DetPhs(di,wi,fi);
end

DataPhs = DataPhs.';

% Absorb the old SD.SrcAmp/DetPhs coefficients into my new SrcPhs and
% DetPhs coefficients (not needed in DataPhs because genBornData()
% already included the old SD.SrcPhs SD.DetPhs values).

SrcPhs = SrcPhs .* SD.SrcAmp;
DetPhs = DetPhs .* SD.DetAmp;

return;
