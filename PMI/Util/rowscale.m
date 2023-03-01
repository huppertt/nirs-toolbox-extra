%ROWSCALE       Scale the rows of a matrix or array.
%
%   As = rowscale(A, c)
%
%   As          The row scaled array.
%
%   A           The input array.
%
%   c           The scaling coefficients.  This have same number of rows as
%               A and one column
%
%   Calls: none.
%
%   Bugs: only handles arrays upto 3 dimension.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/02/07 16:30:57 $
%
%  $Revision: 1.2 $
%
%  $Log: rowscale.m,v $
%  Revision 1.2  2001/02/07 16:30:57  dboas
%  removed the 3rd dimension, as this was screwing up our sparse matrix
%  calculation.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function As = rowscale(A, c)

[m n] = size(A);
[cm] = size(c);

if cm ~= m
    error('A and c must have the same number of rows')
end

As = zeros(size(A));
for ir = 1:m
  As(ir, :) = A(ir,:) * c(ir);
end
