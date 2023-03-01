% TD_GRADGF  Gradient of the time domain Green's function for 
%            DPDW (see Patterson, 1989).  
%
% phi = TD_GradGF(r1, r2, t1, t2, D, idxLambda);
%
% Compute gradient of fluence at [r2, t2] given a delta-function
% source at [r1,t1].  r1, r2, t1, t2 can all be vectors.  The size of
% r1,r2 and t1,t2 must agree.  Sizes of r1,t1 and r2,t2 need not agree.
%
% Uses a single D value for all points.

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

function[grad_phi] = TD_GradGF(r1, r2, t1, t2, Medium, iwvl)

c0 = 2.99792458e10;

V = c0 ./ Medium.idxRefr;

D = V(iwvl) ./  (3 * (Medium.Muspo(iwvl) + Medium.Muao(iwvl)));

% Make the two Nx3 matrices r1 and r2 self-consistant (if possible)

if     (size(r1,2) ~= 3) | (size(r2,2) ~= 3)
   error('r1 and r2 must be Nx3 matrices');
elseif (size(r1,1) == 1) & (size(r2,1) ~= 1)
%   disp('Expanding r2');
   r1 = ones(size(r2,1),1) * r1;
elseif (size(r2,1) == 1) & (size(r1,1) ~= 1)
%   disp('Expanding r2');
   r2 = ones(size(r1,1),1) * r2;
elseif (size(r1,1) ~= size(r2,1))
   error('Vectors r1 and r2 have different sizes');
end

% The vector interval, now an Nx3 array

rvec = r2 - r1;

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
   warning('Vectors t1 and t2 have different sizes');
   keyboard;
end

% dt is automatically a Nx1 row vector

dt = t2 - t1;

% TD_GradGF() returns a (ntime x npos) matrix.

phi = TD_GF(r1, r2, t1, t2, Medium, iwvl);

tpos = find(dt >  0);
tneg = find(dt <= 0);

% Gradient is easily computed in terms of the original fluence
%    \grad{\phi} = \frac{-\vec{r}}{2 D t}\phi

if (~isempty(tpos))
   % Must be (ntime x npos) to match TD_GradGF()
   gph(tpos,:) = -1 ./ (2 * D * dt(tpos) * ones(1,size(rvec,1)));
end

if (~isempty(tneg))
   gph(tneg,:) = 0;    % Best I can manage.  Really, it's infinite.
end

tmp_gphi = phi .* gph;          % gradient sans \vec{r} -> ntime by npos

% I can't matrix multiply N-dimensional matrices, do it the hard way.

for k = 1:3
   % Put back \vec{r} component-by-component
   grad_phi(:,:,k) = tmp_gphi .* (rvec(:,k) * ones(1,size(tmp_gphi,1)))';
end

return;
