% GENMUJACOBIAN   Generate the Jacobian for unknown background
%                 Optical Properties
%
% J = genMuJacobian(SD, Medium, MeasList, Method, dMu, MuFlag);
%   
%   SD,Medium,MeasList - PMI fields
%
%   Method - This specifies the method used to generate the Jacobian
%               'Born' - Solves the least squares difference.
%               'Rytov' - Solves the least squares log-difference.
%
%   MuFlag - indicates whether to generate Jacobian for unknown
%      scattering, absorption, or both.
%
%      MuFlag(1) indicates whether to include absorption
%      MuFlag(2) indicates whether to include scattering
%
%   dMu - this is the increment in scattering and absorption for
%      the finite difference estimate of the Jacobian
%
%      dMu(1,:) are the increment for absorption
%      dMu(2,:) are the increment for scattering
%
% Returns:
%   J - the Jacobian given the specified flags (wavelength changes fastest)

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

function [J] = genMuJacobian(SD, Medium, MeasList, Method, dMu, MuFlag)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

%%
%% Error Check
%%
if ~strcmpi(Method,'Born') & ~strcmpi(Method,'Rytov')
   error('ERROR: genSDJacobian must be a method of Born or Rytov!');
end

%%
%%  Loop over measurements
%%

nMeas = size(MeasList,1);
nWvl  = length(SD.Lambda);

J = zeros(nMeas, 2*nWvl);

idxRow = 1;

Phi0 = DPDWHelmholtz(SD, Medium, MeasList);

Mua0  = Medium.Muao;
Musp0 = Medium.Muspo;

%%
%%  Calc Mua Jacobian
%%

if (MuFlag(1))
   %% Mua increment
   
   Medium.Muao = Medium.Muao + dMu(1,:);
   
   Phi1 = DPDWHelmholtz(SD, Medium, MeasList);

   Medium.Muao = Mua0;
   
   B = Phi1 - Phi0;
   
   if strcmp(Method, 'Rytov')
      B = B ./ Phi0;
   end

   Ba = spalloc(nMeas, nWvl, nMeas);
   
   for iWvl = 1:nWvl
      ml = find(MeasList(:,4) == iWvl);
      
      if (~isempty(ml))
	 Ba(ml,iWvl) = B(ml) / dMu(1, iWvl);
      end
   end
else
   Ba = spalloc(nMeas, nWvl, 1);
end
   
%%
%%  Calc Musp Jacobian
%%

if MuFlag(2)
   %% Musp increment

   Medium.Muspo = Medium.Muspo + dMu(2,:);
   
   Phi1 = DPDWHelmholtz(SD, Medium, MeasList);
   
   Medium.Muspo = Musp0;
   
   B = Phi1 - Phi0;

   if strcmp(Method, 'Rytov')
      B = B ./ Phi0;
   end	
   
   Bs = spalloc(nMeas, nWvl, nMeas);
   
   for iWvl = 1:nWvl
      ml = find(MeasList(:,4) == iWvl);
      
      if (~isempty(ml))
	 Bs(ml,iWvl) = B(ml) / dMu(2, iWvl);
      end
   end
else
   Bs = spalloc(nMeas, nWvl, 1);
end

J = [ Ba Bs ];

return;
