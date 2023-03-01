%calcDelMuSp     Calculate the mu_sp perturbation due to the anomalies.
%
%   delMuSp = calcDelMuSp(Objects, CompVol, BkgdMuSp, idxLambda)
%
%   delMuSp      The scattering function perturbation at each voxel.
%
%   Objects     The cell array of PMI Object data structures describing
%               each anomaly.
%
%   CompVol     The PMI Computational Volume structure defining the voxels
%               at which del mu_a will be calculated.
%
%   BkgdMuSp     The background mu_sp parameter for the wavelength.
%
%   idxLambda   OPTIONAL: The wavelength index used to select the
%               absorption parameter for the objects.
%
%
%   Calls: GenSphere, GenBlock
%
%   Bugs: overlapping objects currently produce the sum of their absorption
%   parameters, currently only handles uniform voxelations.

% Copyright (C) 2002, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function delMuSp = calcDelMuSp(Objects, CompVol, BkgdMuSp, idxLambda)

if nargin < 4
    idxLambda = 1;
end
if strcmp(CompVol.Type, 'uniform') == 0
    error('This routine only handles uniform voxelations');
end
%%
%%  pre-allocations
%%
delMuSp = zeros(length(CompVol.Y), length(CompVol.X), length(CompVol.Z));

%%
%%  Loop over the cell array of objects
%%
nObjects = length(Objects);

for iObj = 1:nObjects
    switch lower(Objects{iObj}.Type)
     case 'sphere'
      delMuSp = delMuSp + GenSphere(CompVol, Objects{iObj}.Pos, ...
                                  Objects{iObj}.Radius, ...
                                  Objects{iObj}.Mu_sp(idxLambda));
     case 'block'
      delMuSp = delMuSp + GenBlock(CompVol, Objects{iObj}.Pos, ...
                                  Objects{iObj}.Dims, ...
                                  Objects{iObj}.Mu_sp(idxLambda));
     otherwise
      error(['Unknown object: ' Objects{iObj}.Type]);
    end
end

%%
%%  Subtract the background absoprtion
%%
idxNZ = delMuSp > 0;
delMuSp(idxNZ)  =  delMuSp(idxNZ) - BkgdMuSp(idxLambda);
