% TD_2PTSB  "Exact" solution of Helmholtz equation for a slab.
%
%   Phi0 = TD_2ptSB(SD, Medium, MeasList, [Debug])
%
%   Phi0     The incident response at each detector from each source in a
%            column vector.  The same combination pattern as the columns of A.
%
%  SD, Medium   PMI Structures
%  MeasList     Measurement List
%
%   Debug       OPTIONAL: Print out debugging info.
%
%   TD_2ptSB computes the exact solution of the time-dependant
%   Helmholtz Equation for a spatial uniform D.  A slab boudary
%   condition is implemented using image charges.

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

function[Phi0] = TD_2ptSB(SD, Medium, MeasList, Debug)

if (~exist('MeasList','var')) | (isempty(MeasList))
   MeasList = SD.MeasList;
end

if ~exist('Debug','var')
   Debug = 0;
end

stepping       = 51;
dt0            = 50e-12;
trapezoid_rule = 1;
fixed_interval = 1;

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

nMeas = size(MeasList,1);
Phi0  = zeros(nMeas, 1);

% Individual zBnd calculations are slow, minimize the number needed
%
% Errors in zBnd are the largest source of errors in the integrand
%  (at about the 0.1% level, compared with 1% due to time gates).

if (exist('mexzbnd','file') == 3)
   % Use the good routine, if availabl
   zBnd = mexzbnd(Medium.idxRefr, Medium.Muspo);
else
   zBnd = calcExtBnd(Medium.idxRefr, Medium.Muspo);
end

%%
%% Extract the slab geometry.  
%%

Thickness = Medium.Slab_Thickness;

Z0 = min(0, Thickness);
Z1 = max(0, Thickness);

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
   
   Src = SD.SrcPos(MeasList(iMeas, 1),:);
   Det = SD.DetPos(MeasList(iMeas, 2),:);
   iWv = MeasList(iMeas,4);
   
   v   = V(iWv);
   mus = Medium.Muspo(iWv);
   mua = Medium.Muao(iWv);
   
   %%  Move the effective source position 1 mean free path into the medium

   Src = moveSrcSlab(Src, Medium.Slab_Thickness, mus);

   % Move detectors?  I don't think so, but ....

   Det = moveSrcSlab(Det, Medium.Slab_Thickness, mus);
   
   %%
   %%  Compute the incident response at the detector.
   %%
   
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
	 DT = T2 / (stepping - 1);
      end

      if (T1 > 1e-6)
	 warning('Delay time exceeds 1us, unit conversion problem?');
	 keyboard;
      end
   
      if (T2 > 1e-6)
	 warning('Gate width exceeds 1us, unit conversion problem?');
	 keyboard;
      end
   
      phi = 0;

      Tdet = T1 + [0:DT:T2]';
   else
      % Return instantaneous value if T2 is non-positivie
      Tdet = T1;
   end
   
   T0   = ones(size(Tdet)) * T0;

   % Initialize to source position
   zi1 = Src(1,3);
   zi2 = Src(1,3);
   chg = 1;               % Source has charge +1

   % Get field due to source charge

   phi_src = TD_GF(Src, Det, T0, Tdet, Medium, iWv);
   phi_img = 0;

   bflag   = 0;
   iSrc    = Src;

   if (phi_src <= 0) 
      if (Debug)
	 disp(sprintf('Phi0(%d), phi_src is zero',iMeas));
      end
      
      % If the source charge is zero (e.g. large separations or
      % times), then no number of image charges will make it non-zero.
      
      Phi0(iMeas) = 0;
      continue;
   end
   
   % for all image charges, sum up the field

   for loop = 2:50           
      % Max 50 image charges
      % Get the images of the images at the opposite boundary
      chg = -chg;
      
      zii1 = getImageCharge(zi2, Z0 - sign(Thickness)*zBnd(iWv));
      zii2 = getImageCharge(zi1, Z1 + sign(Thickness)*zBnd(iWv));

      % Fluence due to image charges
      
      iSrc(1,3) = zii1;
      phi_img1 = TD_GF(iSrc, Det, T0, Tdet, Medium, iWv);
      
      iSrc(1,3) = zii2;
      phi_img2 = TD_GF(iSrc, Det, T0, Tdet, Medium, iWv);
      
      % Check for convergence, set flag

      if all(phi_src==0) & all(phi_img1==0) & all(phi_img2==0)
	 % For negative times, fluence is always zero
	 bflag = 1;
      else
	 if (any(phi_src+phi_img)==0)
	    warning('Divide by zero');
	    keyboard;
	 end
	 
	 dphi = mean(abs(phi_img1 + phi_img2)) ./ mean(abs(phi_src + phi_img));
	 
	 if (dphi < 1e-5)
	    bflag = 1;
	 end
	 
	 clear dphi
      end
      
      % Add in contribution for this pair of image charges

      phi_img = phi_img + chg * (phi_img1 + phi_img2);

      % Stop if converged, else copy back and iterate
      
      if (bflag)
	 break;
      else
	 zi1 = zii1;
	 zi2 = zii2;
	 clear zii1 zii2 phi_img1 phi_img2;
      end
   end

   if (~bflag)
      warning('Fluence failed to converge');
      keyboard;
   end

   % This still needs to be back-ported to ZB and NB!!!!
   phitmp = V(iWv) * (phi_src + phi_img);

   if (T2 > 0)
      if (trapezoid_rule | fixed_interval)
	 % Integrate with trapezoid rule
	 phitmp = sum(phitmp(1:end-1) + phitmp(2:end))*DT/2;
      else
	 % Evaulte integral using Simpson's rule.
	 phitmp = (sum(phitmp(1:2:end-1) + phitmp(3:2:end))/2/3 + ...
		   sum(phitmp(2:2:end-1)) * 2/3) * 2*DT;
      end
   else
      % If T2 <= 0, return the instantaneous value, but that's exactly what
      % we've calculated so there's nothing more to do.
   end
   
   sa = SD.SrcAmp(MeasList(iMeas, 1), iWv);
   da = SD.DetAmp(MeasList(iMeas, 2), iWv);

   Phi0(iMeas) = phitmp * sa * da;

   clear phi_src phi_img phitmp
end

if Debug
   fprintf('\n');
end

return;
