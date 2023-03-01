% TIK  Tikhonov regularized inverse
%
% X = tik(A, Y, alpha, [VD], [VI]);
%
% Y     - residues
% A     - forward matrix
% alpha - regularization parameter, multiplies VI (may be a vector)
% VD    - OPTIONAL: variance of the data (or empty list for identity)
% VI    - OPTIONAL: variance of the image (or empty list for identity)
%
% X     - reconstructed image

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

function[X] = tik(A, Y, alpha, VD, VI);

% Sanity check initial arguments

if (size(Y,1) == 1 & size(Y,2) ~= 1)
   % ' or .' ???
   Y = Y.';
end

if (size(A,1) ~= length(Y))
   error('A and Y have different sizes');
end

% Fill in remaining arguments as needed

if (~exist('VD','var') | isempty(VD))
   % Create default identity matrix
   VD = speye(size(A,1));
elseif (length(VD) == size(A,1) | length(VD) == 1)
   % Expand into diagonal matrix
   VD = spdiags(VD(:), 0, size(A,1), size(A,1));
else
   if (any(size(VD) ~= size(A,1)))
      error('Data variance and forward matrix sizes disagree');
   end
end

if (~exist('VI','var') | isempty(VI))
   % Create default identity matrix
   VI = speye(size(A,2));
elseif (length(VI) == size(A,2) | length(VI) == 1)
   % Expand into diagonal matrix
   VI = spdiags(VI(:), 0, size(A,2), size(A,2));
else
   if (any(size(VI) ~= size(A,2)))
      error('Image covariance and forward matrix sizes disagree');
   end
end

% Invert the data using Tikhonov regularization

% X = (A^\dagger A + R)^{-1} (A^\dagger Y)

VI = inv(VI);                    % Need inverse covariance, not covariance
VD = inv(VD);

for n = 1:length(alpha)
   X(:,n) = (A' * VD * A + alpha(n) * VI) \ (A' * VD * Y);
end

return;
