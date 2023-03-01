% GENSDJACOBIAN   Generate the Jacobian for unknown SD Amplitude/Phase
%
% J = genSDJacobian(SD, Medium, MeasList, data, Method);
%   
%   Method - This specifies the method used to generate the Jacobian
%               'Born'  - Solves the least squares difference.
%               'Rytov' - Solves the least squares log-difference.
%
%   SD, Medium - PMI data structures
%
%   data     - This is the measured fluence for each source-detector pair.
%              data is required for creating the sdAmp Jacobian using
%              the Born method, for the Rytov method it is unused.
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

function[J] = genSDJacobian(SD, Medium, MeasList, Phi, Method)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

%%
%% Error Check
%%
if ~strcmpi(Method,'Born') & ~strcmpi(Method,'Rytov')
   error('genSDJacobian must be a method of Born or Rytov!');
end

if (strcmpi(Method,'Born') & (length(Phi) ~= size(MeasList,1)))
   error('Phi and MeasList must have equal lengths');
end

%%
%%  Determine various numbers 
%%
nWvl  = max(1,length(SD.Lambda));
nFrq  = max(1,length(SD.ModFreq));
nSrc  = size(SD.SrcPos,1);
nDet  = size(SD.DetPos,1);
nMeas = size(MeasList,1);

%% How much space do we need?

nCC = nWvl * nFrq * (nSrc + nDet);

%%
%% Add the SD Amplitudes as unknowns in the Jacobian
%%

J = zeros(nMeas,nCC);

for iMeas = 1:nMeas
   iSrc = MeasList(iMeas,1);
   iDet = MeasList(iMeas,2);
   iFrx = MeasList(iMeas,3);
   iWvx = MeasList(iMeas,4);
      
   iFrq = max(1,iFrx);
   iWvl = max(1,iWvx);
   
   iJS = (iFrq-1)*nSrc*nWvl + (iWvl-1)*nSrc + iSrc;
   iJD = (iFrq-1)*nDet*nWvl + (iWvl-1)*nDet + iDet + nSrc*nWvl*nFrq;
   
   if (strcmpi(Method,'Rytov'))
      J(iMeas, iJS) = 1;
   else
      J(iMeas, iJS) = Phi(iMeas) / SD.SrcAmp(iSrc, iWvl, iFrq);
   end

   if (strcmpi(Method,'Rytov'))
      J(iMeas, iJD) = 1;
   else
      J(iMeas, iJD) = Phi(iMeas) / SD.DetAmp(iDet, iWvl, iFrq);
   end
end

J = sparse(J);

return;

