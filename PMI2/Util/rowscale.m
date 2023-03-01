% ROWSCALE  Scale the rows or columns of a matrix or array.
%
% As = rowscale(A, s);
%
%   A           The input array.
%
%   c           The scaling coefficients.  If the length of c is equal
%               to the number of rows of A, then scale the rows of A by
%               s, otherwise scale its columns by s.
%
% Bugs: only handles arrays upto 3 dimension.

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

function[A] = rowscale(A, s)

[n1,n2] = size(A);
[nc]  = length(s);

% Matrix multiplies should be much faster than a for loop and should yield
% the same answer in the end.

if     (nc == n1)
   % for iRow = 1:n1
   %    A(iRow, :) = A(iRow,:) * s(iRow);
   % end

   A = spdiags(s, 0, n1, n1) * A;
elseif (nc == n2)
   % for iCol = 1:n2
   %    A(:,iCol) = A(:,iCol) * s(iCol);
   % end
   
   A = A * spdiags(s, 0, n2, n2);
else
   error('A and s have different sizes');
end

return;

