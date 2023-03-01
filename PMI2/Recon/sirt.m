% SIRT   Simultaneous iterative reconstruction technique
%
% X = sirt(A, Y, X0, nIter, W);
%
%   X     The estimate of the x vector.
%   A     The forward matrix.
%   Y     The measured data.
%
%   X0    OPTIONAL: An initial guess, if not supplied then 0 will be
%         used.  You can use [] as a placeholder if needed.
%
%   nIter OPTIONAL: The maximum number of iterations to compute
%         (default 10 * number of rows).  Use [] as a placeholder.
%
%   W     OPTIONAL: The relaxation parameter, defining how far to "step"
%         on each iteration.  A value of 1 causes each step to reach
%         the hyperplane of the current orthogonal projection.  Less than
%         1 cause the step to fall short of the hyperplane.  The default
%         value is 1.  Use [] as a placeholder.

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


function[X] = sirt(A, phi, X0, nIter, w)

Debug = 0;

[nmeas, nvox] = size(A);

if (size(phi,1) ~= nmeas)
   error('Size of A and phi do not agree.  Transposed data or matrix?');
end

%%
%%  Default arguments.  Allow empty matrices [] to be used as
%%   place-holders.  Same defaults as if they were never specified in
%%   the first place.
%%

if (~exist('w','var') | isempty(w))
   w = 1;
end

if (~exist('nIter','var') | isempty(nIter))
   nIter = 10;
end

if (~exist('X0','var') | isempty(X0))
   X0 = zeros(nvox, 1);
end

nEst = length(nIter);

X = zeros(nvox, nEst);

%%
%%  Precompute the row norm
%%

if (Debug)
   disp('Computing row norms');
end

rownorm = zeros(nmeas,1);

for i = 1:nmeas
   % Was
   %    rownorm(i) = A(i,:) * A(i,:)';
   % but this will be much faster for sparse matrices
   
%   Arow = full(A(i,:));
   Arow = A(i,:);
   rownorm(i) = Arow * Arow';
end

zlst = find(rownorm==0);

if (~isempty(zlst))
   rownorm(zlst) = 9.0e50;  %% deal with zero signal rows - DAB 99-10-12
end

clear zlst

%%
%%  Loop over the number of iterations requested
%%

if (Debug)
   disp('Reconstructing Image');
end

for j = 1:nEst
   %%
   %%  Copy the estimate from the previous iteration parameter to the
   %%  estimate for the new x(:,j)iteration parameter.
   %%

   if (j == 1)
      X(:, j) = X0;
      N       = nIter(1);
   else
      disp(sprintf('Starting estimate %d', j));
      
      X(:,j) = X(:,j-1);
      N      = nIter(j) - nIter(j-1);
   end

   %%
   %%  Iterate over each row
   %%
   
   xnew = zeros(nvox, 1);
   
   for i = 1:N
      RelResid = (phi - A * X(:,j)) ./ rownorm;
      xnew = full(RelResid' * A)' / nmeas;
      
      X(:,j) = X(:,j) + w * xnew;

      %%  Add positivity constraint some day
   end
end

return;

