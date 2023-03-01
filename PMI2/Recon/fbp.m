% FBP   Filtered back-projection
%
% X = fbp(A, Y, alpha, [VD], [VI]);
%
% Y     - residues
% A     - forward matrix
% alpha - regularization parameter, multiplies CI (may be a vector)
% VD    - OPTIONAL: variance of the data (or empty list for identity)
% VI    - OPTIONAL: variance of the image (or empty list for identity)
%
% X     - reconstructed image
%
% VD and VI are _ignored_ until I can figure out the proper scalings.

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

function[X] = fbp(A, Y, alpha, VD, VI);

% Sanity check initial arguments

if (size(Y,1) == 1 & size(Y,2) ~= 1)
   % ' or .' ???
   Y = Y.';
end

if (size(A,1) ~= length(Y))
   error('A and Y have different sizes');
end

% Invert the data using filtered back-projection

% X = A^\dagger (A A^\dagger + R)^{-1} Y

for n = 1:length(alpha)
   X(:,n) = A' * ((A * A' + alpha(n) * speye(size(A,1))) \ Y);
end

return;
