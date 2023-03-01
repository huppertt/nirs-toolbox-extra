% TDGF  Time domain Green's function for DPDW (see Patterson, 1989)
%
% phi = TD_GF(r1, r2, t1, t2, Medium, idxLambda);
%
% Compute fluence at [r2, t2] given a delta-function source at
%  [r1,t1].  r1, r2, t1, t2 can all be vectors.  The size of r1,r2
%  and t1,t2 must agree.  Sizes of r1,t1 and r2,t2 need not agree.
%
% Uses a single D value for all measurements

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

function[phi] = TD_GF(r1, r2, t1, t2, Medium, idxLambda);

c0 = 2.99792458e10;
V  = c0 ./ Medium.idxRefr(idxLambda);

% Note that this is different from Patterson's definition of D
D  = V / (3 * (Medium.Muspo(idxLambda) + Medium.Muao(idxLambda)));

% Make the two Nx3 matrices r1 and r2 self-consistant (if possible)

if     (size(r1,2) ~= 3) | (size(r2,2) ~= 3)
   warning('r1 and r2 must be Nx3 matrices');
   keyboard
elseif (size(r1,1) == 1) & (size(r2,1) ~= 1)
%   disp('Expanding r2');
   r1 = ones(size(r2,1),1) * r1;
elseif (size(r2,1) == 1) & (size(r1,1) ~= 1)
%   disp('Expanding r2');
   r2 = ones(size(r1,1),1) * r2;
elseif (size(r1,1) ~= size(r2,1))
   error('Vectors r1 and r2 have different sizes');
end

% Compute the square-distance between points.  Make dR^2 a row vector

dR2 = sum((r2 - r1).^2, 2)';

% Time can also be a vector

if     (size(t1,2) ~= 1) | (size(t2,2) ~= 1)
   warning('t1 and t2 must be Nx1 matrices');
   keyboard;
elseif (size(t1,1) == 1) & (size(t2,1) ~= 1)
%   disp('Expanding t2');
   t1 = ones(size(t2,1),1) * t1;
elseif (size(t2,1) == 1) & (size(t1,1) ~= 1)
%   disp('Expanding t2');
   t2 = ones(size(t1,1),1) * t2;
elseif (size(t1,1) ~= size(t2,1))
   error('Vectors t1 and t2 have different sizes');
end

% dt is automatically a row vector

dt  = t2 - t1;

% Enforce causality even if the diffusion equation doesn't.  Also
%  handle the exp(-1/0) -> 0 case.

lc = find(dt > 0);            % Causal times

A0     = zeros(size(dt));
A0(lc) = 1 ./ sqrt((4*pi * D * dt(lc)).^3);
A0     = A0 * ones(size(dR2));

% 1e6 to kill non-causal times
B0     = 1e6*ones(size(dt));
B0(lc) = 1./(4 * D * dt(lc));

C0     = V * Medium.Muao(idxLambda) * dt;
C0     = C0 * ones(size(dR2));

phi    = A0 .* exp(-B0 * dR2 - C0);

return;
