% GENIDXJACOBIAN   Generate the Jacobian for unknown background
%                  index of refraction
%
% J = genIdxJacobian(SD, Medium, MeasList, Method);
%   
%   SD,Medium,MeasList - PMI fields
%
%   Method - This specifies the method used to generate the Jacobian
%               'Born' - Solves the least squares difference.
%               'Rytov' - Solves the least squares log-difference.
%
% Returns:
%   J - the Jacobian given the specified flags (wavelength changes fastest)

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

function [J] = genIdxJacobian(SD, Medium, MeasList, Method)

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

idx0  = Medium.idxRefr;
dN    = 0.1;  % Increment

%%
%%  Calculate Jacobian
%%

Medium.idxRefr = Medium.idxRefr + dN;
   
Phi1 = DPDWHelmholtz(SD, Medium, MeasList);

Medium.idxRefr = idx0;
   
B = Phi1 - Phi0;
   
if strcmp(Method, 'Rytov')
   B = B ./ Phi0;
end

J = spalloc(nMeas, nWvl, nMeas);
   
for iWvl = 1:nWvl
   ml = find(MeasList(:,4) == iWvl);
      
   if (~isempty(ml))
      J(ml,iWvl) = B(ml) / dN;
   end
end

return;
