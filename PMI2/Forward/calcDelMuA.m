% CALCDELMUA   Calculate the mu_a perturbation due to the anomalies.
%
% delMuA = calcDelMuA(Medium, Objects, idxLambda)
%
%   Medium      The PMI structure ...
%
%   Objects     The cell array of PMI Object data structures describing
%               each absorption anomaly [OPTIONAL]
%
%   idxLambda   The wavelength index used to select the absorption
%               parameter for the objects.
%
% Returns:
%   delMuA      The absorption function perturbation at each voxel.
%
% Bugs: overlapping objects currently produce the sum of their absorption
%   parameters

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

function [delMuA] = calcDelMuA(Medium, Objects, idxLambda)

if (~exist('Objects') | isempty(Objects))
   Objects = Medium.Object;
end

%% I would prefere to use a Measurement list, but what's needed
%%  is indeed the index into the Mu_a values, which depend only
%%  on the CompVol.

if (~exist('idxLambda','var'))
   idxLambda = 1;
end

delMuA = [];

%%
%%  Loop over the cell array of objects
%%
if (~isempty(Objects))
   nObjects = length(Objects);

   for iObj = 1:nObjects
      if (isfield(Objects{iObj}, 'idxRefr'))
         warning([ 'Object{' num2str(iObj) '}.idxRefr ignored' ]);
      end

      switch lower(Objects{iObj}.Type)
	 case 'sphere'
	    dMuA = genSphere(Medium, Objects{iObj}.Pos, ...
			     Objects{iObj}.Radius, ...
			     Objects{iObj}.Mua(idxLambda) ...
					- Medium.Muao(idxLambda));
	 case 'block'
	    dMuA = genBlock(Medium, Objects{iObj}.Pos, ...
			    Objects{iObj}.Dims, ...
			    Objects{iObj}.Mua(idxLambda) ...
				       - Medium.Muao(idxLambda));
	 case 'image'
	    dMuA = Objects{iObj}.Mua(:,idxLambda) ...
		                       - Medium.Muao(idxLambda);
	 otherwise
	    error(['Unknown object: ' Objects{iObj}.Type]);
      end
      
      if (isempty(delMuA))
	 delMuA = dMuA;
      else
	 delMuA = delMuA + dMuA;
      end
   end
end

return;

