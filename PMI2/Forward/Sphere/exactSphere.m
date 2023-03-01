% EXACTSPHERE  Analytical solution for DPDW scattering off of a sphere
%               in an infinite medium
%
% phiS = exactSphere(SD, Medium, MeasList, Object);
%
% SD, Medium, MeasList - usual PMI structures
% Object               - Object structure for this perturbation
%

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

function[phiS] = exactSphere(SD, Medium, MeasList, Object);

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (isTD(SD, MeasList))
   error('exactSphere() only supports CW/RF measurements');
end

if (~isstruct(Object) | ~strcmpi(Object.Type,'sphere'))
   error('Object must be a sphere');
end

if (size(Object.Pos,1) == 3 & size(Object.Pos,2) == 1)
   % Make Object.Pos match SD.SrcPos and SD.DetPos.
   Object.Pos = Object.Pos';
end

iFrq = unique(MeasList(:,3));
iWvl = unique(MeasList(:,4));

if (length(iFrq) ~= 1)
   error('MeasList cannot have more than one frequency');
end

if (length(iWvl) ~= 1)
   error('MeasList cannot have more than one wavelength');
end

% Highest order spherical harmonic (l) to sum over
lmax = 20;

% Geometrical constants

Rsph = Object.Radius;

rsrc = SD.SrcPos(MeasList(:,1),:);
rdet = SD.DetPos(MeasList(:,2),:);
opos = ones(size(MeasList,1),1) * Object.Pos;

rso = sqrt(sum((rsrc - opos).^2,2));
rod = sqrt(sum((opos - rdet).^2,2));
% rsd = sqrt(sum((rsrc - rdet).^2,2));

% Cosine of the angle between src-obj and obj-det 
%  (computed as the dot product of the unit vectors)

r1 = (opos - rsrc) ./ (rso * ones(1,3));
r2 = (opos - rsrc) ./ (rod * ones(1,3));

c_theta = sum(r1 .* r2, 2);

clear rsrc rdet opos r1 r2;

% Optical properties

vo  = 2.99792458e10 / Medium.idxRefr(iWvl);
w   = 2 * pi * (SD.ModFreq(iFrq) * 1e6);

Do  = vo / (3 * (Medium.Muspo(iWvl) + Medium.Muao(iWvl)));
Ko  = sqrt((-vo * Medium.Muao(iWvl) + i * w) / Do);

% Object.idxRefr is not supported elsewhere, make it an optional field
if (isfield(Object,'idxRefr'))
   vi = 2.99792458e10 / Object.idxRefr(iWvl);
else
   vi = vo;
end

Di  = vi / (3 * (Object.Musp(iWvl) + Object.Mua(iWvl)));
Ki  = sqrt((-vi * Object.Mua(iWvl) + i * w) / Di);
	 
% Compute solution

x = Ko * Rsph;
y = Ki * Rsph;
  
phiS = zeros(size(MeasList,1),1);

for l = [0:lmax];
   % Compute all Y_{l,m} for this l
   SHs = SphHarm(l, -1 * ones(size(c_theta)), 0);
   SHd = SphHarm(l,                  c_theta, 0);
   
   % l-dependant part
   A = -i * Ko * bes_h1(l,Ko*rso) * ((Do * x *  bes_jd(l,x) *  bes_j(l,y) -...
				      Di * y *   bes_j(l,x) * bes_jd(l,y))/...
				     (Do * x * bes_h1d(l,x) *  bes_j(l,y) -...
				      Di * y *  bes_h1(l,x) * bes_jd(l,y)));
   
   % Sum over the m-dependant parts
   
   phi_l = SHs(1) * SHd(1);
    
   for m = 1:l
      phi_l = phi_l + conj(SHs(m+1)) * SHd(m+1);
      phi_l = phi_l + SHs(m+1) * conj(SHd(m+1));
   end

   % Update the scattered field
   
   phiS = phiS + A .* bes_h1(l, Ko * rod) * phi_l;
end
  
%  phi_3pt = (vo/Do) * (exp(j*Ko*rso+j*Ko*rod)/(4*pi*rso*4*pi*rod)) * 
%                         (4/3*pi*a^3*vo/Do);

phiS = vo/Do * phiS;

return;
 
