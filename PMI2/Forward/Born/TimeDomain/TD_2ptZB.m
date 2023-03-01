% TD_2PTZB  "Exact" solution of Helmholtz equation for a semi-infinte volume
%
%   Phi0 = TD_2ptZB(SD, Medium, MeasList, Debug)
%
%   Phi0     The incident response at each detector from each source in a
%            column vector.  The same combination pattern as the columns of A.
%
%  SD, Medium   PMI Structures
%  MeasList     Measurement List
%
%   Debug       OPTIONAL: Print out debugging info.
%
%   TD_2ptZB computes the exact solution of the time-dependant
%   Helmholtz Equation for a spatial uniform D in a semi-infinte
%   medium.  The extrapolated boundary condition is met by using image
%   charges.

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

function[Phi0] = TD_2ptZB(SD, Medium, MeasList, Debug)

if (~exist('MeasList','var')) | (isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var'))
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

if nargin < 4
   Debug = 0;
end

nMeas = size(MeasList,1);
Phi0  = zeros(nMeas, 1);

if (exist('mexzbnd','file') == 3)
   zBnd = mexzbnd(Medium.idxRefr, Medium.Muspo);
else
   zBnd = calcExtBnd(Medium.idxRefr, Medium.Muspo);
end


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
   iWvl = MeasList(iMeas, 4);
   
   %%  Move the effective source position 1 mean free path into the medium

   if (~isfield(Medium,'Slab_Thickness') | Medium.Slab_Thickness > 0)
      rSrc(3) = max(0,rSrc(3)) + 1./Medium.Muspo(iWvl);
      rDet(3) = max(0,rDet(3)) + 1./Medium.Muspo(iWvl);
   else
      rSrc(3) = min(0,rSrc(3)) - 1./Medium.Muspo(iWvl);
      rDet(3) = min(0,rDet(3)) - 1./Medium.Muspo(iWvl);
   end

   % Locate image charge
   
   rImg = rSrc;
   rImg(3) = getImageCharge(rSrc(3), -sign(Medium.Slab_Thickness)*zBnd(iWvl));

   %%
   %%  Compute the incident response at the detector.
   %%

   % Source pulse (a delta-function) defines t=0
   
   % Source pulse (a delta-function) defines t=0

   T0 = 0 + SD.SrcOffset(MeasList(iMeas,1),MeasList(iMeas,4));
   T1 = SD.TimeDelay(MeasList(iMeas, 6)) + ...
	SD.DetOffset(MeasList(iMeas,2),MeasList(iMeas,4));
   T2 = SD.TimeGateWidth(MeasList(iMeas, 7));

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

      % Integrate time-gate using trapezoid rule
   
      phi_o = 0;
      phi_i = 0;

      Tdet = T1 + [0:DT:T2]';
      T0   = ones(size(Tdet)) * T0;
   
      if (isempty(Tdet))
	 % Trivial
      
	 Phi0(iMeas) = 0;
	 continue;
      end
   
      % TD_GF is vectorized
      phi_o = V(iWvl) * TD_GF(rSrc, rDet, T0, Tdet, Medium, iWvl);
      phi_i = V(iWvl) * TD_GF(rImg, rDet, T0, Tdet, Medium, iWvl);

      if (trapezoid_rule | fixed_interval)
	 % Integrate with trapezoid rule
	 
	 phitmp = sum(phi_o(1:end-1) + phi_o(2:end))*DT/2;
	 phi_o  = phitmp;
   
	 phitmp = sum(phi_i(1:end-1) + phi_i(2:end))*DT/2;
	 phi_i  = phitmp;
      else
	 % Evaulte integral using Simpson's rule.
	 phitmp = (sum(phi_o(1:2:end-1) + phi_o(3:2:end))/2/3 + ...
		   sum(phi_o(2:2:end-1)) * 2/3) * 2*DT;
	 phi_o = phitmp;

	 phitmp = (sum(phi_i(1:2:end-1) + phi_i(3:2:end))/2/3 + ...
		   sum(phi_i(2:2:end-1)) * 2/3) * 2*DT;
	 phi_i = phitmp;
      end
   else
      % Return instantaneous value if width is non-positive
      
      phi_o = V(iWvl) * TD_GF(rSrc, rDet, T0, T1, Medium, iWvl);
      phi_i = V(iWvl) * TD_GF(rImg, rDet, T0, T1, Medium, iWvl);
   end
   
   sa = SD.SrcAmp(MeasList(iMeas, 1), iWvl);
   da = SD.DetAmp(MeasList(iMeas, 2), iWvl);

   Phi0(iMeas) = (phi_o - phi_i) * sa * da;
end

if Debug
   fprintf('\n');
end

return;


