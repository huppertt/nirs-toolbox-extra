% GENJACOBIANFROM2PTS  Compute the Mu Jacobian for FullBorn matrices
%
% J = genJacobianFrom2pts(SD, Medium, MeasList, ...
%                         Phi2pt, PhiTot, Method, OptProp);
%
%   SD,Medium,MeasList - PMI fields
%
%   Phi2pt, PhiTot     - returned by genBornData() using 'FullBorn'
%
%   OptProp            - optical property flags
%      OptProp(1) indicates whether to include absorption
%      OptProp(2) indicates whether to include scattering
%
%   Method             - method to use to calculate J, 'Born' or 'Rytov'
%
% Returns:
%   J - the Jacobian given the specified flags.

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

function[J] = genJacobianFrom2pts(SD, Medium, MeasList, ...
				  Phi2pt, PhiTot, Method, OptProp, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

% Check the arguments

switch lower(Method)
   case 'born'
      PhiTot = ones(size(MeasList,1), 1);
   case 'rytov'
      if (~exist('PhiTot','var') | isempty(PhiTot))
	 error('PhiTot not defined for Rytov method (required)');
      end
   otherwise
      error([ 'Unsupported method ' Method ]);
end

% Recreate the mapping used to generate Phi2pt and PhiTot.  This will
% likely fail if the measurement list has changed since Phi2pt and
% PhiTot were generated, but I don't see any way around that.

FB = MLtoFB(SD, Medium, MeasList, OptProp);

% Compute the matrix

% Phi2pt is (nVoxs+nDets)*nOpts by (nSrc+nDets)*nOpts
% PhiTot is already in measurement list order.
%
% What is this SUPPOSED to be?  Should be the change in each
% measurement due to an infinitesimal change in the optical properties
% of each voxel.  If so, I really need the Src->Voxel terms too, don't I?

J = zeros(size(MeasList,1), size(Phi2pt,1));

for iMeas = 1:size(MeasList,1)
   if (Debug)
      disp(num2str(iMeas));
   end
   
   iOpt = FB.optR(iMeas);   % Optical property set
   iSrc = FB.srcR(iMeas);   % Src index
   iDet = FB.detR(iMeas);   % Det index

   phi_sv = Phi2pt(:,iSrc);
   phi_dv = Phi2pt(:,FB.nSrcs + iDet);
   
   J(iMeas,:) = (phi_sv .* phi_dv).' / PhiTot(iMeas);
   clear phi_sv phi_dv;
end

% Repack from nvox*nwvl stack to standard Born ordering (block-matrix form)

J = repackBornFwdMat(SD, MeasList, J, OptProp, Debug);

return;
