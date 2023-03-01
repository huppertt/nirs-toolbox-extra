% GENPOSJACOBIAN  Generate the Jacobian for unknown SD Positions
%
% J = genPosJacobian(SD, Medium, MeasList, Method, posFlags);
%   
%   Method - This specifies the method used to generate the Jacobian
%               'Born'  - Solves the least squares difference.
%               'Rytov' - Solves the least squares log-difference.
%
%   SD, Medium - PMI data structures
%
%   posFlags - is a flag with 6 elements. 1-Yes, 0-No.
%         sdPos(1)   indicates whether to include the X positions of the srcs
%         sdPos(2)   ... Y position of sources ...
%         sdPos(3)   ... Z position of sources ...
%         sdPos(4)   ... X position of detectors ...
%         sdPos(5)   ... Y position of detectors ...
%         sdPos(6)   ... Z position of detectors ...
%
% Returns:
%   J - the Jacobian given the specified flags.

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

function[J] = genPosJacobian(SD, Medium, MeasList, Method, sdPos)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

%%
%% Error Check
%%
if ~strcmpi(Method,'Born') & ~strcmpi(Method,'Rytov')
   error('genSDJacobian must be a method of Born or Rytov!');
end

if (sdPos(3) ~= 0 | sdPos(6) ~= 0)
   if (strcmpi(Medium.Geometry,'slab'))
      % Changing Z involves changing the slab thickness
      if (Debug)
	 warning('Slab thickness implicitly changed in slab geometries');
      end
   end
   
   if (strcmpi(Medium.Geometry,'semi-infinite') | ...
       strcmpi(Medium.Geometry,'semi'))
      % I can't do Z-motion in a semi-infinite geometry
      warning('Unable to move Z positions in semi-infinite geometries');
      sdPos(3) = 0;
      sdPos(6) = 0;
   end
   
   % Infinite doesn't need a warning
end

%%
%%  Determine various numbers 
%%
nWvl  = length(SD.Lambda);
nSrc  = size(SD.SrcPos,1);
nDet  = size(SD.DetPos,1);
nMeas = size(MeasList,1);

%%
%%  Add positional uncertainty to matrix
%%

J = zeros(nMeas, 3*(nSrc + nDet));

Phi0 = DPDWHelmholtz(SD, Medium, MeasList);

%% Source uncertainty in X
if sdPos(1)
   B = doPosCorr(SD, Medium, MeasList, Phi0, Method, 1, 'X');
   J(:,0*nSrc + [1:nSrc]) = B;
end

%% Source uncertainty in Y
if sdPos(2)
   B = doPosCorr(SD, Medium, MeasList, Phi0, Method, 1, 'Y');
   J(:,1*nSrc + [1:nSrc]) = B;
end

%% Source uncertainty in Z
if sdPos(3)
   B = doPosCorr(SD, Medium, MeasList, Phi0, Method, 1, 'Z');
   J(:,2*nSrc + [1:nSrc]) = B;
end

%% Detector uncertainty in X
if sdPos(4)
   B = doPosCorr(SD, Medium, MeasList, Phi0, Method, 0, 'X');
   J(:,3*nSrc + 0*nDet + [1:nDet]) = B;
end

%% Detector uncertainty in Y
if sdPos(5)
   B = doPosCorr(SD, Medium, MeasList, Phi0, Method, 0, 'Y');
   J(:,3*nSrc + 1*nDet + [1:nDet]) = B;
end

%% Detector uncertainty in Z
if sdPos(6)
   B = doPosCorr(SD, Medium, MeasList, Phi0, Method, 0, 'Z');
   J(:,3*nSrc + 2*nDet + [1:nDet]) = B;
end

if (any(sdPos==0))
   J = sparse(J);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[B] = doPosCorr(SD, Medium, MeasList, Phi0, Method, isSrc, gradDir)

%% dR is the stepping used to calculate the gradients

dR = 0.01;

%% System size

nMeas = size(MeasList,1);
nSrc  = size(SD.SrcPos,1);
nDet  = size(SD.DetPos,1);
nWvl  = length(SD.Lambda);

%% Output size

if (isSrc)
   B = zeros(nMeas, nSrc);
else
   B = zeros(nMeas, nDet);
end

% Displacements to use

if (isSrc)
  dRsrc = dR;
  dRdet =  0;
else
  dRsrc =  0;
  dRdet = dR;
end

% Direction to displace

if     ((gradDir >= 'X') & (gradDir <= 'Z'))
  iDir = gradDir - 'X' + 1;
elseif ((gradDir >= 'x') & (gradDir <= 'z'))
  iDir = gradDir - 'x' + 1;
else
  warning([ 'Unknow direction ' gradDir ]);
  keyboard;
end

% Displace the positions

if (iDir < 3 | ~strcmpi(Medium.Geometry, 'slab'))
   SD.SrcPos(:,iDir) = SD.SrcPos(:,iDir) + dRsrc;
   SD.DetPos(:,iDir) = SD.DetPos(:,iDir) + dRdet;
else
   % Adjust slab thickness
   
   if (dRsrc == 0)
      SD.DetPos(:,3) = SD.DetPos(:,3) + dRdet;
      Medium.Slab_Thickness = Medium.Slab_Thickness + dRdet;
   else
      % Keep srcs at 0, move detectors to compensate
      
      SD.DetPos(:,3) = SD.DetPos(:,3) - dRsrc;
      Medium.Slab_Thickness = Medium.Slab_Thickness - dRsrc;
   end
end

% Calculate the perturbed fluence

Phi2 = DPDWHelmholtz(SD, Medium, MeasList);

% Compute Jacobian elements (pull outside loop over measurements so improve
% efficiency of vector operation)

dPhi = (Phi2 - Phi0) / dR;

if (strcmpi(Method, 'Rytov'))
   dPhi = dPhi ./ Phi0;
end

% Repack as matrix

for iMeas = 1:nMeas;
   iSrc = MeasList(iMeas,1);
   iDet = MeasList(iMeas,2);
   iWvl = MeasList(iMeas,4);

   if (isSrc)
      iB = iSrc;
   else
      iB = iDet;
   end
   
   B(iMeas, iB) = dPhi(iMeas);
end

return;

