% GENSPHEREDATA   Generate field scattered off a sphere in an infinite
%                  medium using the analytical solution.
%
% Phi = genSphereData(SD, Medium, MeasList);

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

function[Phi] = genSphereData(SD, Medium, MeasList, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

if (~isfield(Medium,'Object'))
   error('Medium.Object doesn''t exist - no perturbing objects!');
end

if (isTD(SD, MeasList))
   error('Analytical sphere solution only supports CW/RF measurements');
end

%%
%%  Generate measured data
%%

if (~strcmpi(Medium.Geometry, 'infinite') & ...
    ~strcmpi(Medium.Geometry, 'inf')      & ...
    ~strcmpi(Medium.Geometry, 'inft'))
   error('Analytical spherical solution only valid in infinite media');
end
  
%%
%%  Check the object(s) to make sure that they are spheres
%%

nObj = length(Medium.Object);

if (nObj > 1)
   warning([ 'Analytical solution for multiple spheres only valid in' ...
	     ' first Born approximation' ]);
end

for k = 1:nObj
   if ~strcmpi(Medium.Object{k}.Type, 'sphere')
      error([ 'Medium.Object{' num2str(k) '} does not define a sphere' ]);
   end
end
  
%% Calculate Perturbations

PhiS = zeros(size(MeasList,1), 1);

vfrq = unique(MeasList(:,3));
vwvl = unique(MeasList(:,4));

for iWvl = 1:length(vwvl)
   if (Debug)
      disp(sprintf('Wavelength -> %d nm', SD.Lambda(iWvl)));
   end
   
   for iFrq = 1:length(vfrq)
      if (Debug)
	 disp(sprintf('Modulation Frequency -> %d MHz', SD.ModFreq(iFrq)));
      end
      
      ml = find(MeasList(:,3) == vfrq(iFrq) & MeasList(:,4) == vwvl(iWvl));
      
      if (~isempty(ml))
	 for iObj = 1:nObj
            if (Debug)
               disp([ 'Ojbect ' num2str(iObj) ]);
            end

	    % Compute exact solution for each sphere
	    phiObj = exactSphere(SD, Medium, ...
				 MeasList(ml,:), Medium.Object{iObj});

	    % Apply SD amplitudes
	    sAmp = SD.SrcAmp(MeasList(ml,1),iWvl);
	    dAmp = SD.DetAmp(MeasList(ml,2),iWvl);
	    
	    % Sum scattered field(s)
	    PhiS(ml) = PhiS(ml) + phiObj .* sAmp .* dAmp;
	    
	    clear sAmp dAmp phiObj;
	 end
      end
   end
end

% Calculate unperturbed field
Phi0 = DPDWHelmholtz(SD,Medium,MeasList);

% Return total fluence
Phi = Phi0 + PhiS;

return;
