% FD_GF  Frequency-domain two-point Green's function
%
% Phi = FD_GF(rSrc, rVox, D, K, v);
%
%   rSrc - 1x3 vector, source position
%   rVox - Nx3 vectors, voxel positions
%   D    - Diffusion coefficient
%   K    - Wave-vector
%   v    - Speed of light
%
% Returns:
%   Phi  - The incident response at each voxel

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

function [phi] = FD_GF(rSrc, rVox, D, K, v)

% rSrc and rVox must be position vectors

[ns1, ns2] = size(rSrc);
[nv1, nv2] = size(rVox);

if (ns2 ~= 3)
   error('rSrc is not an Nx3 position vector');
end

if (nv2 ~= 3)
   error('rVox is not an Nx3 position vector');
end

%% Source and "detector" must have the same size

if (ns1 == 1 & nv1 ~= 1)
   rSrc = ones(nv1,1) * rSrc;
end

if (nv1 == 1 & ns1 ~= 1)
   rVox = ones(ns1,1) * rVox;
end

%%
%%  Compute the incident response at the "detector".
%%

rSV = sqrt(sum((rVox - rSrc).^2,2));

if (any(rSV <= 1e-25))
   warning('Zero-separation not allowed');
   keyboard;
end

% Green's function for photon density (freqency-domain)

phi = exp(-K * rSV) ./ (4 * pi * D * rSV);

return;
