% ART  Algebraic Reconstruction Technique
%
% X = art(A, Y, X0, nIter, W, iMeas)
%
%   X     The estimate of the reconstructed data
%   A     The forward matrix.
%   Y     The measured data.
%
%   X0    OPTIONAL: An initial guess, if not supplied then 0 will be
%         used.  Use [] for a placeholder.
%
%   nIter OPTIONAL: The maximum number of iterations to compute.
%         If nIter is a vector an estimate is returned for each 
%         element in nIter.  The number iterations must be increasing.  
%         Use [] for a placeholder.
%
%   W     OPTIONAL: The relaxation parameter, defining how far to "step"
%         on each iteration.  A value of 1 causes each step to reach
%         the hyperplane of the current orthogonal projection.  Less than
%         1 cause the step to fall short of the hyperplane.  The default
%         value is 1.  Use [] for a placeholder.
%
%   iMeas OPTIONAL: The initial row to project onto.  Useful for
%         examining convergence performance.  The default value is 1.

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


function[X] = art(A, phi, X0, nIter, w, iMeas)

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
   nIter = 10*nvox;
end

if (~exist('X0','var') | isempty(X0))
   X0 = zeros(nvox, 1);
end

if (~exist('iMeas','var') | isempty(iMeas))
   iMeas = 1;
end

nEst = length(nIter);

X = zeros(nvox, nEst);

%%
%%  Precompute the row norm
%%

rownorm = zeros(nvox,1);

for i = 1:nmeas
   rownorm(i) = A(i,:) * A(i,:)';
end

zlst = find(rownorm==0);

if (~isempty(zlst))
   %% deal with zero signal rows - DAB 99-10-12
   rownorm(zlst) = 9.0e50;  
end

clear zlst

%% Loop over measurements in random order or the update tends to become
%% spatially correlated.

pmeas = randperm(nmeas);

%%
%%  Loop over the number of iterations requested
%%

for j = 1:nEst
   %%
   %%  Copy the estimate from the previous iteration parameter to the
   %%  estimate for the new x(:,j)iteration parameter.
   %%

   if (j == 1)
      X(:,j) = X0;
      N      = nIter(1);
   else
      disp(sprintf('Starting estimate %d',j));

      X(:,j) = X(:,j-1);
      N      = nIter(j) - nIter(j-1);
   end

   %%
   %%  Iterate over each row 
   %%

   for i = 1:N
      jMeas = pmeas(iMeas);
      
      RelResid = (phi(jMeas) - A(jMeas,:) * X(:,j)) ./ rownorm(jMeas);
      xnew     = (RelResid * A(jMeas,:))';

      X(:,j)   = X(:,j) + w * xnew;

      iMeas = mod((iMeas-1)+1, nmeas) + 1;

      %%  Add positivity constraint
   end
end

return;

