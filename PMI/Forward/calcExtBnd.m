%calcExtBnd     Calculate the extrapolated boundary distance.
%
%   dBnd = calcExtBnd(idxRefr, mu_sp)
%
%   dBnd        The distance to the extrapolated boundary.
%
%   idxRefr     The index of refraction for the diffuse medium.
%
%   mu_sp       The reduced scattering coefficient for the diffuse medium.
%
%
%   calcExtBnd returns the distance to the extrapolated boundary from a
%   air-diffuse medium interface.  This is calculated using the approach
%   described in:
%
%    Boundary conditions for the diffusion equation in radiative transfer
%    Haskell et. al.,  J. Opt. Soc. Am - A   Vol. 11, No. 10,  Oct 1994
%    pg 2727 - 2741
%
%
%   Calls: none.
%
%   Bugs: none known.

% Copyright (C) 2002, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: calcExtBnd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.1  1999/12/03 13:55:15  dboas
%  This routine now does the exact calculation for the extrapolated
%  boundary condition.  This calculation is taken from the
%  Haskel paper.
%
%  Revision 3.0  1999/06/17 17:39:56  rjg
%  Initial Revision for PMI 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dBnd = calcExtBnd(idxRefr, mu_sp)


oc = asin(1/idxRefr);

ostep = pi / 2000;

o = 0:ostep:oc;

cosop = (1-idxRefr^2 * sin(o).^2).^0.5;
coso = cos(o);
r_fres = 0.5 * ( (idxRefr*cosop-coso)./(idxRefr*cosop+coso) ).^2;
r_fres = r_fres + 0.5 * ( (idxRefr*coso-cosop)./(idxRefr*coso+cosop) ).^2;

r_fres(ceil(oc/ostep):1000) = 1;

o = 0:ostep:ostep*(length(r_fres)-1);
coso = cos(o);

r_phi_int = 2 * sin(o) .* coso .* r_fres;
r_phi = sum(r_phi_int) / 1000 * pi/2;

r_j_int = 3 * sin(o) .* coso.^2 .* r_fres;
r_j = sum(r_j_int) / 1000 * pi/2;

Reff = (r_phi + r_j) / (2 - r_phi + r_j);
dBnd = 2/3 * (1 + Reff) / (1 - Reff) / mu_sp;
