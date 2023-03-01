% TD_2PTNB  "Exact" solution of Helmholtz equation (infinite volume)
%
%   Phi0 = TD_2ptNB(SD, Medium, MeasList, Debug)
%
%   Phi0     The incident response at each detector from each source in a
%            column vector.  The same combination pattern as the columns of A.
%
%  SD, Medium   PMI Structures
%  MeasList     Measurement list
%
%   Debug       OPTIONAL: Print out debugging info.
%
%   TD_2ptNB computes the exact solution of the time-dependant
%   Helmholtz Equation for a spatial uniform D in an infinite
%   volume.

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

function[Phi0] = TD_2ptNB(SD, Medium, MeasList, Debug)

if (~exist('MeasList','var')) | (isempty(MeasList))
   MeasList = SD.MeasList;
end

if ~exist('Debug','var')
   Debug = 0;
end

stepping       = 51;
dt0            = 10e-12;
trapezoid_rule = 0;
fixed_interval = 0;
      
%
% Local Speed of light
%

V = 2.99792458e10 ./ Medium.idxRefr;

% Per-optode temporal shifts

if (~isfield(SD,'SrcOffset') | isempty(SD.SrcOffset))
   SD.SrcOffset = zeros(size(SD.SrcPos,1),length(SD.Lambda));
end

if (~isfield(SD,'DetOffset') | isempty(SD.DetOffset))
   SD.DetOffset = zeros(size(SD.DetPos,1),length(SD.Lambda));
end

% Discretization used in calculating the convlution
DT = min(SD.TimeGateWidth)/25;

nMeas = size(MeasList,1);
Phi0  = zeros(nMeas, 1);

%%
%%  Loop over each source-detector combination
%%

if Debug
   fprintf('Meas #: ');
end

for iMeas = 1:nMeas
   if Debug
      fprintf('%d  ', iMeas);
   end

   rSrc = SD.SrcPos(MeasList(iMeas, 1),:);
   rDet = SD.DetPos(MeasList(iMeas, 2),:);
   iwvl = MeasList(iMeas, 4);
   
   %%
   %%  Compute the incident response at the detector.
   %%

   % Source pulse (a delta-function) defines t=0

   T0 = 0 + SD.SrcOffset(MeasList(iMeas,1),MeasList(iMeas,4));
   T1 = SD.TimeDelay(MeasList(iMeas, 6)) + ...
	SD.DetOffset(MeasList(iMeas,2),MeasList(iMeas,4));
   T2 = SD.TimeGateWidth( MeasList(iMeas, 7));

   % Shift times so that the pulse once again arrives at T = 0

   T1 = T1 - T0;
   T0 = 0;

   if (T2 > 0)
      % Discretization used in calculating the convlution

      if (fixed_interval)
	 DT = dt0;
      else
	 DT = min(SD.TimeGateWidth)/(stepping - 1);
      end

      % Integrate over time gate using trapezoid rule
      
      phi  = 0;
      Tdet = T1 + [0:DT:T2]';
      T0   = ones(size(Tdet)) * T0;
      
      % TD_GF is vectorized
      phi = V(iwvl) * TD_GF(rSrc, rDet, T0, Tdet, Medium, iwvl);

      if (trapezoid_rule | fixed_interval)
	 % Integrate with trapezoid rule
	 phitmp = sum(phi(1:end-1) + phi(2:end))*DT/2;
      else
	 phitmp = (sum(phi(1:2:end-1) + phi(3:2:end))/2/3 + ...
		   sum(phi(2:2:end-1)) * 2 / 3) * 2*DT;
      end
   else
      % Return instantaneous value
      
      phitmp = V(iwvl) * TD_GF(rSrc, rDet, T0, T1, Medium, iwvl);
   end
   
   sa = SD.SrcAmp(MeasList(iMeas, 1), iwvl);
   da = SD.DetAmp(MeasList(iMeas, 2), iwvl);

   Phi0(iMeas) = phitmp * sa * da;
   
   clear phi phitmp sa da
end

if Debug
   fprintf('\n');
end

return;


