% TD_3PTBORN1ZB   1st Born approx. of a perturbed Helmholtz equation.
%
%   TD_3PTBORN1NB computes the forward weighting matrix associated
%   with the Born-1 approximation to a spatial varying k in a 
%   slab geometry.
%
% [Phi0,A] = TD_3ptBorn1SB(SD, Medium, MeasList, OptProp, Debug);
%
%  SD,Medium,MeasList  PMI structures
%  OptProp             Flags, OptProp(1)->abs, OptProp(2)->scat
%  Debug               OPTIONAL: Print out debugging info.
%
% Returns:
%   A        The forward matrix
%   Phi0     The unperturbed fluence
%
% Bugs:
%   This code hasn't been tested nearly as well as the (faster and more 
%   useful) MEX versions.  While it gives more or less correct answers,
%   it may still have bugs lurking in the odd corners.

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

function[Phi0, A] = TD_3ptBorn1SB(SD, Medium, MeasList, OptProp, Debug)

if (~exist('MeasList','var')) | (isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var'))
   Debug = 0;
end

V = 2.99792458e10 ./ Medium.idxRefr;
D  = V ./ (3 * (Medium.Muspo + Medium.Muao));

nMeas = size(MeasList,1);

if (length(OptProp(1) < 2))
   % Pad out options vector
   OptProp(end+1:2) = 0;
end

if (OptProp(1) ~= 0)
  calc_mua = 1;
else
  calc_mua = 0;
end

if (OptProp(2) ~= 0)
  calc_musp = 1;
else
  calc_musp = 0;
end

%%
%%  Create the sampling volume
%%

[Xm Ym Zm dV] = sampleVolume(Medium.CompVol);

rVox = [ Xm, Ym, Zm ];
nvox = length(Xm);
clear Xm Ym Zm;

%%
%% Matrix weighting
%%

dmua = -1 * ones(size(V));
dmus =  D ./ V;

%% Allocate storage space for the matrix early on

if calc_mua & calc_musp
  A = zeros(nMeas, 2*nvox);
else
  A = zeros(nMeas,  nvox);
end

%% Calculate extrapolated boundaries just once

if (exist('mexzbnd','file') == 3)
   % Use the more accurate compiled form if available.
   zBnd = mexzbnd(Medium.idxRefr, Medium.Muspo);
else
   zBnd = calcExtBnd(Medium.idxRefr, Medium.Muspo);
end

% Calculate the incident fluence, use existing 2-pt routines

if (Debug)
   disp('Calculating incident fluence');
end

Phi0 = TD2pt(SD, Medium, MeasList, Debug);

% If OptProp(:) == 0, then the matrix is trivial

if (calc_mua == 0 & calc_musp == 0)
   A = [];        % Nothing more to do
   return;
end

%%
%%  Loop over each source detector combination.  
%%  Memory for A has already been allocated.
%%

stepping = 25;

if Debug
   fprintf('Meas #: ');
end

for iMeas = 1:nMeas
   if Debug
      fprintf('%d  ', iMeas);
   end

   %%  Compute the incident response on every voxel from each measurement 
   %%  (source detector pair).  The source is represented by a point source.

   iSrc = MeasList(iMeas, 1);
   iDet = MeasList(iMeas, 2);
   iWvl = MeasList(iMeas, 4);

   sA = SD.SrcAmp(iSrc, iWvl);
   dA = SD.DetAmp(iDet, iWvl);
   
   if (isempty(iSrc) | isempty(iDet) | isempty(iWvl))
      error('Index not found');
   end

   rSrc = SD.SrcPos(iSrc,:);
   rDet = SD.DetPos(iDet,:);

   % rVox alread defined
   
   %% Convolve sources, detectors, and gates
   
   T0 = 0;
   T1 = SD.TimeDelay(     MeasList(iMeas, 6) );
   T2 = SD.TimeGateWidth( MeasList(iMeas, 7) );

   DT   = T2 / (stepping-1);
   Tdet = T1 + [0:DT:T2]';

   Aa = zeros(length(Tdet),size(rVox,1));
   As = zeros(length(Tdet),size(rVox,1));

   % Use the analytical solutions for the imaging kernels

   if (calc_mua)
      Aa = GetPhiA(rSrc, rDet, rVox, T0, Tdet, Medium, V, D, iWvl, zBnd);
   end      

   if (calc_musp)
      As = GetPhiS(rSrc, rDet, rVox, T0, Tdet, Medium, V, D, iWvl, zBnd);
   end

   % Evaulte integral dTdet using trapezoid rule.
   
   if (calc_mua)
      tmpA = sum(Aa(1:end-1,:) + Aa(2:end,:)) * DT/2;
   end
   if (calc_musp)
      tmpS = sum(As(1:end-1,:) + As(2:end,:)) * DT/2;
   end
   
   if (calc_mua & calc_musp)
      A(iMeas,1:nvox)           = tmpA * dmua(iWvl) * sA * dA;
      A(iMeas+nMeas,nvox+1:end) = tmpS * dmus(iWvl) * sA * dA;
   elseif (calc_mua)	        
      A(iMeas,1:nvox)           = tmpA * dmua(iWvl) * sA * dA;
   elseif (calc_musp)	        	      
      A(iMeas,1:nvox)           = tmpS * dmus(iWvl) * sA * dA;
   end
   
   clear tmpA tmpS
end

clear Aa As iSrc iDet

% dV because perturbation is proportional to the size of the voxel,

A = A * dV;

if Debug
   fprintf('\n');
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute analytical solution for the absorption imaging kernel for
%  a vector of times and positions

function[A] = GetPhiA(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl, zBnd)

if (size(rs,1) ~= 1 | size(rd,1) ~= 1)
   warning('rs and rd must be 1x3 vectors');
end

% Move "sources" one scattering length into the volume

rs = moveSrcSlab(rs, Medium.Slab_Thickness, Medium.Muspo(iWvl));
rd = moveSrcSlab(rd, Medium.Slab_Thickness, Medium.Muspo(iWvl));

% Volume bounds (upper and lower limits to slab)

Z1 = min(0, Medium.Slab_Thickness) - zBnd(iWvl);
Z2 = max(0, Medium.Slab_Thickness) + zBnd(iWvl);

% Pre-compute the image charges

zs(1,1:2) = rs(3);

for k = ([1:50] + 1)
   zs(k,1) = getImageCharge(zs(k-1,2), Z1);
   zs(k,2) = getImageCharge(zs(k-1,1), Z2);
end

zd(1,1:2) = rd(3);

for k = ([1:50] + 1)
   zd(k,1) = getImageCharge(zd(k-1,2), Z1);
   zd(k,2) = getImageCharge(zd(k-1,1), Z2);
end

% Allocate/initialize work vectors

rsi = rs;  rsi(3) = zs(1,1);
rdi = rd;  rdi(3) = zd(1,1);

% Start with infinite-medium solution

A0 = GetAa_inft(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl);

if (all(A0(:) == 0))
   % If A0 is zero, no number of image charges will help

   A = A0;
   return;
end

for loop = 1:50
   A1 = zeros(size(A0));
   
   % Not efficient Matlab, but at least it's clearly correct

   for il = -loop:loop
      if (il >= 0)
	 rsi(3) = zs( il+1,1);
      else
	 rsi(3) = zs(-il+1,2);
      end
	 
      for jl = -loop:loop
	 if (abs(il) ~= loop & abs(jl) ~= loop)
	    % Add on only the new terms, terms with il~=loop AND jl~=loop
            % have already been computed.
	    
	    continue;
	 end

	 if (jl >= 0)
	    rdi(3) = zd( jl+1,1);
	 else
	    rdi(3) = zd(-jl+1,2);
	 end
	 
	 % Sign of this image charge
	 
	 if (mod(abs(il+jl),2) == 0)
	    chg = 1;
	 else
	    chg = -1;
	 end
	 
	 % Contribution of this pair of "sources"
	 
	 A1 = A1 + ...
	      chg * GetAa_inft(rsi, rdi, rv, T0, Tdet, Medium, V, D, iWvl);
      end
   end
      
   % Stop if we've converged sufficiently
      
   if (sum( A1(:)./(A0(:) + A1(:)) ) < 1e-5)
      A = A0 + A1;
      return;
   end
   
   % If A is still zero after the update, no number of charges will help

   if (all(A1(:) == 0.0) & loop > 1)
      A = A0 + A1;
      return;
   end
   
   % Copy back updated result and iterate
   
   A0 = A0 + A1;
end

if (loop >= 50)
   error('Failed to converge');
end

A = A0;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute analytical solution for the scattering imaging kernel for
%  a vector of times and positions

function[A] = GetPhiS(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl, zBnd)

if (size(rs,1) ~= 1 | size(rd,1) ~= 1)
   warning('rs and rd must be 1x3 vectors');
end

% Move "sources" one scattering length into the volume

rs = moveSrcSlab(rs, Medium.Slab_Thickness, Medium.Muspo(iWvl));
rd = moveSrcSlab(rd, Medium.Slab_Thickness, Medium.Muspo(iWvl));

% Volume bounds (upper and lower limits to slab)

Z1 = min(0, Medium.Slab_Thickness) - zBnd(iWvl);
Z2 = max(0, Medium.Slab_Thickness) + zBnd(iWvl);

% Pre-compute the image charges

zs(1,1:2) = rs(3);

for k = ([1:50] + 1)
   zs(k,1) = getImageCharge(zs(k-1,2), Z1);
   zs(k,2) = getImageCharge(zs(k-1,1), Z2);
end

zd(1,1:2) = rd(3);

for k = ([1:50] + 1)
   zd(k,1) = getImageCharge(zd(k-1,2), Z1);
   zd(k,2) = getImageCharge(zd(k-1,1), Z2);
end

% Allocate/initialize work vectors

rsi = rs;  rsi(3) = zs(1,1);
rdi = rd;  rdi(3) = zd(1,1);

% Start with infinite-medium solution

A0 = GetAs_inft(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl);

if (all(A0(:) == 0))
   % If A0 is zero, no number of image charges will help

   A = A0;
   return;
end

for loop = 1:50
   A1 = zeros(size(A0));
   
   % Not efficient Matlab, but at least it's clearly correct

   for il = -loop:loop
      if (il >= 0)
	 rsi(3) = zs( il+1,1);
      else
	 rsi(3) = zs(-il+1,2);
      end
	 
      for jl = -loop:loop
	 if (abs(il) ~= loop & abs(jl) ~= loop)
	    % Add on only the new terms, terms with il~=loop AND jl~=loop
            % have already been computed.
	    
	    continue;
	 end

	 if (jl >= 0)
	    rdi(3) = zd( jl+1,1);
	 else
	    rdi(3) = zd(-jl+1,2);
	 end
	 
	 % Sign of this image charge
	 
	 if (mod(abs(il+jl),2) == 0)
	    chg = 1;
	 else
	    chg = -1;
	 end
	 
	 % Contribution of this pair of "sources"
	 
	 A1 = A1 + ...
	      chg * GetAs_inft(rsi, rdi, rv, T0, Tdet, Medium, V, D, iWvl);
      end
   end
      
   % Stop if we've converged sufficiently
      
   if (sum( A1(:)./(A0(:) + A1(:)) ) < 1e-5)
      A = A0 + A1;
      return;
   end
   
   % If A is still zero after the update, no number of charges will help

   if (all(A1(:) == 0.0) & loop > 1)
      A = A0 + A1;
      return;
   end
   
   % Copy back updated result and iterate
   
   A0 = A0 + A1;
end

if (loop >= 50)
   error('Failed to converge');
end

A = A0;

return;

% ----------------------------------------------------------------------
% Infinite medium solution - absorbing kernel

function[A] = GetAa_inft(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl);

% Get the vectors I need

r0 = zeros(size(rv,1),3);
r1 = (rv - ones(size(rv,1),1)*rs);
r2 = (rv - ones(size(rv,1),1)*rd);
   
% Magnitudes of r1 and r2

mr1 = sqrt(sum(r1.^2,2));
mr2 = sqrt(sum(r2.^2,2));
   
% Vector with correct magnitude and random orientation

rTmp      = zeros(length(mr1),3);
rTmp(:,3) = mr1 + mr2;

% TD_GF vectorized in both time and space

A = V(iWvl) * TD_GF(r0, rTmp, T0, Tdet, Medium, iWvl);
	 
% Apply geometric scale factors to get sensitivity and return the value

A = A * V(iWvl) / (4 * pi * D(iWvl));
A = A .* ((1./mr1 + 1./mr2) * ones(1,size(A,1))).';

return;

% ----------------------------------------------------------------------
% Infinite medium solution - scattering kernel

function[A] = GetAs_inft(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl);

dt = Tdet - T0;
ml = find(dt > 0);

% Get the vectors I need

r0 = zeros(size(rv,1),3);
r1 = (rv - ones(size(rv,1),1)*rs);
r2 = (rv - ones(size(rv,1),1)*rd);
   
% Magnitudes of r1 and r2

mr1 = sqrt(sum(r1.^2,2));
mr2 = sqrt(sum(r2.^2,2));
   
% Fake vector with correct magnitude and random orientation

rTmp      = zeros(length(mr1),3);
rTmp(:,3) = mr1 + mr2;

% Start with the absorption kernel

A = GetAa_inft(rs, rd, rv, T0, Tdet, Medium, V, D, iWvl);

% Go from abs to scat, avoiding the divide-by-zero cases

s = zeros(size(A));

if (~isempty(ml))
   % Add back the dot product of the unit vectors
   A = A .* (ones(size(A,1),1) * (sum(r1 .* r2, 2) ./ (mr1 .* mr2))');

   % Convert Kabs to Kscat
   s(ml,:) = ones(length(ml),1) * (mr1.^2 - mr1.*mr2 + mr2.^2)' ./ ...
	     (2 * D(iWvl) * dt(ml) * (mr1 .* mr2)');
   s(ml,:) = s(ml,:) + ones(length(ml),1) * (mr1' + mr2').^2 ./ ...
	     ((2 * D(iWvl) * dt(ml)).^2 * ones(1,size(mr1,1)));
   
   A = A .* s;
end

return;

