% TCGLS   Truncated CG algorithm for solving A'Ax = A'b.
%
%   [X,resid] = tcgls(A, Y, X0, nIter, W)
%
%   X     The estimate(s) of the x vector.
%   resid The final residual vector.
%
%   A     The forward matrix.
%   Y     The measured data.
%
%   X0    OPTIONAL: An initial guess, if not supplied then 0 will be
%         used.  You can use [] as a placeholder if needed.
%
%   nIter OPTIONAL: The maximum number of iterations to compute
%         (default is 1 iteration).  Use [] as a placeholder if needed.
%
%   W     OPTIONAL: The relaxation parameter, defining how far to "step"
%         on each iteration.  A value of 1 causes each step to reach
%         the hyperplane of the current orthogonal projection.  Less than
%         1 cause the step to fall short of the hyperplane.  The default
%         value is 1.
%
%   TCGLS computes the truncated conjugate gradient solution on the normal
%   equations without explicitly forming the normal equations.

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

function[X, R] = tcgls(A, Y, X0, nIter, w)

if (~exist('w','var') | isempty(w))
   w = 1.0;
end

if (~exist('nIter','var') | isempty(nIter))
   nIter = 1;
end

if (~exist('X0','var'))
   % isempty(X0) is handled farther down
   X0 = [];
end

%%
%%  Initalizations
%%

[nmeas,nvox] = size(A);
nEst = length(nIter);

if ((nvox > nmeas) & ~issparse(A))
   bproj = A' * Y;
else
   % Transpose gets expensive
   bproj = (Y' * A)';
end

resid = bproj;
d     = bproj;
delta = bproj' * bproj;

%%
%%  Loop over the number of iterations requested
%%

X = zeros(nvox, length(nIter));
R = zeros(nvox, length(nIter));

for j = 1:length(nIter)
   if (j == 1)
      if (isempty(X0))
         X0 = zeros(nvox,1);
	 % else use user-supplied starting point
      end
      
      N = nIter(1);
   else
      X0 = X(:,j-1);
      
      N = nIter(j) - nIter(j-1);
   end

   %%
   %% Compute the iterations between the previous number requested and the
   %% current number requested.
   %%
   
   for i = 1:N
      if ((nvox > nmeas) & ~issparse(A))
	 q = A' * (A * d);
      else
         % Avoid the expensive matrix transpose
	 q = ((A * d)' * A)';
      end
      
      alpha = delta / (d' * q);

      X0 = X0 + w * alpha * d;

      resid = resid - w * alpha * q;
      
      delta_old = delta;
      delta     = resid' * resid;
      beta      = delta / delta_old;

      d = resid + beta * d;
   end

   X(:,j) = X0;
   R(:,j) = resid;

   if (length(nIter) > 1)
      disp(sprintf('Completed %d iterations',nIter(j)));
   end
end

return;
