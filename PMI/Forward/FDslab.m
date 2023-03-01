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

%x0 = 0; xf = 6;		% x-min x-max (cm) surface plane
%y0 = 0; yf = 6;		% y-min y-max (cm) surface plane
%z0 = 0; zf = 3;		% z-min z-max (cm) tissue depth
%h = 0.1;			% step size (cm)
%us1 = 10;			% (mus') Reduced scattering coefficient (cm-1)
%ua1 = 0.05;			% Absorption coefficient (cm-1)

function [A, x, y, z] = FD_slab( x0, xf, y0, yf, z0, zf, h, us1, ua1 )

% Three Dimensional solution to the Diffusion Approximation

E = 1;					% Incident irradiance (W/cm^2)
J = round((xf-x0)/h+1);			% Number of rows in matrix A
K = round((yf-y0)/h+1);			% Number of columns in matrix A
L = round((zf-z0)/h+1);			% Number of whozits in matrix A
midJ = round(J/2);				% middle row
midK = round(K/2);				% middle column

x = x0:h:xf; y = y0:h:yf; z = z0:h:zf;

utr1 = ua1 + us1;		% Total attenuation coefficient (cm-1)
ueff1 = sqrt(3*ua1*utr1);	% Effective attenuation coefficient (cm-1)
D1 = 1/(3*utr1);		% Diffusion coefficient
C1 = -6-ua1*h^2/D1;	% Physical constant for normal tissue
psi0 = 4*E/(1+2*D1*ueff1); % Fluence in cell under laser input

C2 = 1+h/2/D1;			% Surface boundary condition factor
C3 = 0.99;				% Side and bottom boundary condition factor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Coefficient Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Building coefficient matrix...');
f = flops;
G = sparse(J*L*K,J*L*K);

% Define indices for each region
iin 		= GetIndices(2:J-1,	2:K-1,	2:L-1,	J,K,L);
itop 		= GetIndices(1:J,		1:K,		1,			J,K,L);
ileft 	= GetIndices(1:J,		1,			2:L,		J,K,L);
iright 	= GetIndices(1:J,		K,			2:L,		J,K,L);
ifront 	= GetIndices(1,		2:K-1,	2:L,		J,K,L);
iback 	= GetIndices(J,		2:K-1,	2:L,		J,K,L);
ibottom 	= GetIndices(2:J-1,	2:K-1,	L,			J,K,L);

% Define diagonal coefficient matrix
diags = -[0 -1 1 -J J -J*K J*K];
Dmat = zeros(J*K*L,7);
Dmat(iin,1) = C1;
Dmat(iin,2:7) = 1;
Dmat(itop,1) = C2;
Dmat(itop,7) = -1;
Dmat(ileft,1) = -1;
Dmat(ileft,5) = C3;
Dmat(iright,1) = -1;
Dmat(iright,4) = C3;
Dmat(ifront,1) = -1;
Dmat(ifront,3) = C3;
Dmat(iback,1) = -1;
Dmat(iback,2) = C3;
Dmat(ibottom,1) = C2;
Dmat(ibottom,6) = -1;

G(iin(1),iin(1)) = C1;	% Prime matrix with 1 value to overcome spdiags bug
G = spdiags(Dmat,diags,G);
G = G';

disp(sprintf('%s\n','Done'));
Num_Flops = flops-f

% Define H matrix
H = sparse(J*K*L,1);

% Define laser input to surface
j = midJ; k = midK; l = 1;
i = j+(k-1)*J+(l-1)*J*K;
%G(i,:) = sparse(1,J*K*L);
%G(i,i) = 1;
H(i) = psi0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for fluence rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Inverting coefficient matrix...')
f = flops;
s = symmmd(G);			% Creat symmetric minimum degree ordering
spparms('autommd', 0)		% Disable automatic permutations
[ll,uu] = luinc(G(:,s),1e-4);	% Incomplete LU factorization
InitialGuess = ones(J*K*L,1);	% Default initial guess of 0's results in error
psi = bicg(G(:,s),H,1e-6,50,ll,uu,InitialGuess);	% Invert by biconjugate gradient
spparms('autommd',1)				% Enable automatic permutations
psi(s) = psi;						% Reorder solution
Num_Flops = flops-f

A = reshape(psi,J,K,L);

