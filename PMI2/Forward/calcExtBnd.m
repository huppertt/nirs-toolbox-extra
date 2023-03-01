% CALCEXTBND  Calculate the extrapolated boundary distance.
%
% calcExtBnd returns the distance to the extrapolated boundary from a
% air-diffuse medium interface.  This is calculated using the approach
% described in:
%
%  Boundary conditions for the diffusion equation in radiative transfer
%  Haskell et. al.,  J. Opt. Soc. Am - A   Vol. 11, No. 10,  Oct 1994
%  pg 2727--2741.
%
% dBnd = calcExtBnd(idxRefr, mu_sp)
%
%   idxRefr     The index of refraction for the diffuse medium.
%
%   mu_sp       The reduced scattering coefficient for the diffuse medium.
%
% Returns:
%   dBnd        The distance to the extrapolated boundary.

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

function[dBnd] = calcExtBnd(idxRefr, mu_sp)

% The error in the boundary is rougly linear in the number of steps.  This
% used to be 1000, but it turns out errors in the extrapolated boundary are
% probably the single largest source of error in calculating the 2-pt
% functions.

nstep = 10001; % Must be odd for Simpson's rule to work

if (length(idxRefr) ~= length(mu_sp))
   warning('idxRefr and mu_sp have different lengths');
   keyboard;
end

% Initialize to manifestly illegal value
dBnd = nan * ones(size(idxRefr));

% Only calculate for physically meaningful parameters
lidx = find(idxRefr > 0 & mu_sp > 0);

for kidx = 1:length(lidx)
   l = lidx(kidx);
   oc    = asin(1./idxRefr(l));
   ostep = pi / (2*nstep);

   o = 0:ostep:oc;

   cosop = sqrt(1 - (idxRefr(l) * sin(o)).^2);
   coso  = cos(o);

   r_fres =          0.5 * ((idxRefr(l) .* cosop - coso) ./ ...
			    (idxRefr(l) .* cosop + coso)).^2;
   r_fres = r_fres + 0.5 * ((idxRefr(l) .* coso  - cosop) ./ ...
			    (idxRefr(l) .* coso  + cosop)).^2;

   r_fres(ceil(oc/ostep):nstep) = 1;

   o    = 0:ostep:ostep*(length(r_fres)-1);
   coso = cos(o);

   r_phi_int = 2 * sin(o) .* coso .* r_fres;
   r_phi     = sum(r_phi_int(1:end-1) + r_phi_int(2:end)) / (2*nstep) * pi/2;

   r_j_int = 3 * sin(o) .* coso.^2 .* r_fres;
   r_j     = sum(r_j_int(1:end-1) + r_j_int(2:end)) / (2*nstep) * pi/2;

   Reff    = (r_phi + r_j) / (2 - r_phi + r_j);
   dBnd(l) = 2/3 * (1 + Reff) / (1 - Reff) / mu_sp(l);
end

return;
