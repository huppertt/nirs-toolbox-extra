% FD_GRADGF  Frequency-domain two-point Green's function
%
% GradPhi = FD_GradGF(rSrc, rVox, D, K, v);
%
%   rSrc     - 1x3 vector, source position
%   modFreq  - Modulation frequency (in MHz) of this source
%   rVox     - Nx3 vectors, voxel positions
%
% Returns:
%   GradPhi  - The (spatial) gradient of the incident response at each voxel

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

function [GPhi] = FD_GradGF(rSrc, rVox, D, K, v)

[ns1, ns2] = size(rSrc);
[nv1, nv2] = size(rVox);

if (ns1 == 3 & ns2 ~= 3)
   Src = rSrc';
elseif (ns2 == 3)
   Src = rSrc;
else
   warning('rSrc is not a Nx3 position vector');
   keyboard;
end

if (nv1 == 3 & nv2 ~= 3)
   Vox = rVox';
elseif (nv2 == 3)
   Vox = rVox;
else
   warning('rVox is not a Nx3 position vector');
   keyboard;
end

if (size(Src,1) == 1)
   % Source and "detector" must have the same size
   Src = ones(size(Vox,1),1) * Src;
end

%%
%%  Compute the incident response at the "detector".
%%

rvec = Vox - Src;
r    = sqrt(sum(rvec.^2,2));

% Start with original Green's function
Phi0 = FD_GF(rSrc, rVox, D, K, v);

% Analytical deriviative pulls down prefactor and turns it into a vector

% Divide by r here to turn rvec below into a unit vector.
Gtmp = (-K - 1./r) .* Phi0 ./ r;

% This is marginally faster than multiplying by ones and doing a
% single ./ (divison by r=|rvec| was done above)
GPhi = zeros(length(Gtmp),3); % was 8.12 sec
GPhi(:,1) = Gtmp .* rvec(:,1);
GPhi(:,2) = Gtmp .* rvec(:,2);
GPhi(:,3) = Gtmp .* rvec(:,3);

return;
