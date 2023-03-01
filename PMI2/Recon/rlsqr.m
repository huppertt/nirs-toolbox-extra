% RLSQR    Regularized least-squares solution for X
%
%   X = rlsqr(A, Y, lambda, X0, nIter, tol);
%
%   X      The estimate(s) of the x vector.
%
%   A      The forward matrix.
%   Y      The measured data.
%   lambda Regularization param (multiplies identity matrix) (std of noise?)
%
%   nIter  OPTIONAL: The maximum number of iterations to compute.
%          Use [] as a placeholder if needed.
%
%   X0     OPTIONAL: An initial guess, if not supplied then 0 will be
%          used.  You can use [] as a placeholder if needed.
%
%   tol    OPTIONAL: The stopping tolerance.  Default value is 1e-6.
%
% RLSQR calls the matlab routine lsqr() which minimizes the norm
%  of (Y - A*X) with an added amplitude constraint on X.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang,
%                     Jonathan Stott, Greg Boverman
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

function[X] = rlsqr(A, Y, lambda, X0, nIter, tol)

if (~exist('X0','var'))
   % isempty(X0) is handled farther down
   X0 = [];
end

if (~exist('nIter','var') | isempty(nIter))
   nIter = 20;
end

if (~exist('tol','var') | isempty(tol))
   % Generally get good convergence, so keep tol fairly small
   % Is this an absolute or a relative parameter?
   tol = 1e-10;
end

%%
%%  Initalizations
%%

[nmeas,nvox] = size(A);

if (nmeas >= nvox)
   warning('Already have more measurements than voxels - not regularized');
else
   % Extend the matrix

   disp('Extending matrix');

   % speye() gives me classic Tikhonov.  Could just as easily do spatial
   % regulariztaion, laplacian, etc.

   A = [ A; lambda * speye(nvox,nvox) ];

   Y(size(A,1)) = 0;
end

disp('Solving least-squares problem');

if (~isempty(X0))
   X0 = sparse(size(A,2),1);
end

if (~isempty(X0) & any(X0 ~= 0))
   warning('X0 is ignored');
end

[X, tol, niter] = lsqr_greg(A, Y, tol, nIter);

if (niter >= nIter)
   disp(sprintf(...
       'lsqr failed to converge after %d iterations [tol=%g]',niter, tol));
end

return;

% function[solution, tol, num_iter] = lsqr_greg(A, Y, min_tol, max_iter)
% 
% A -> augmented matrix
% Y -> augmented vector 
%
% Greg Boverman
% 7-24-03

function[solution, tol, num_iter] = lsqr_greg(A, Y, min_tol, max_iter)

tol      = 1;
num_iter = 0;

u     = Y;
beta  = norm(u);
u     = (1/beta) * u;

v     = (u' * A)';                          % Avoid taking A'
alpha = norm(v);
v     = (1/alpha) * v;

w = v;

solution = zeros(size(A,2), 1);

phi_overbar = beta;
rho_overbar = alpha;

while ((tol > min_tol) & (num_iter < max_iter))
   u     = A * v - alpha * u;
   beta  = norm(u);
   u     = (1/beta) * u;

   v     = (u' * A)' - beta * v;            % Avoid taking A'
   alpha = norm(v);
   v     = (1/alpha) * v;

   rho = sqrt(rho_overbar.^2 + beta.^2);

   c = rho_overbar / rho;
   s = beta        / rho;

   theta       =  s * alpha;
   rho_overbar = -c * alpha;
   
   phi         =  c * phi_overbar;
   phi_overbar =  s * phi_overbar;

   solution = solution + (phi/rho)*w;

   residnorm = norm(A*solution - Y) / norm(Y);
   w = v - (theta / rho) * w;

   tol = phi_overbar * alpha * abs(c);
   num_iter = num_iter + 1;

   disp(sprintf('%3d: tol=%g', num_iter, tol));
end

return;
