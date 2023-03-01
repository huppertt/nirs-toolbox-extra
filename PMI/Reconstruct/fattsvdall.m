%FATTSVDALL     Truncated SVD solution to a fat matrix problem for all r.
%
%   [xtsvd U S V]= fattsvdall(A, b, U, S, V)
%
%   xtsvd       The estimates of x (n x m), the first column is the estimate
%               at r=1, the last column at r=m.
%
%   U,S,V       The singular value decomposition of A.  If this is supplied
%               on the right hand side it will not be computed.  If only 3
%               input arguments are supplied the economy SVD will be
%               computed.
%
%   A           The system matrix (m x n).  This should not be overdetermined
%               (m <= n).
%
%   b           The measurement vector (m x 1).
%
%
%   FATTSVDALL computed the truncated TSVD solution to an underdetermined
%   system of equations for all truncation parameters.
%
%   Calls: none.
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: fattsvdall.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1999/02/15 19:43:17  rjg
%  Added the ability to handle tall and square problems correctly.
%
%  Revision 1.1  1998/12/22 16:59:18  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xtsvd, U, S, V] = fattsvdall(A, b, U, S, V)

[m n] = size(A);

%%
%%  Compute the economy SVD of the tall problem if necessary
%%
if nargin < 5
    if m < n
        [V S U] = svd(A', 0);
    else
        [U S V] = svd(A, 0);
    end
end
xtsvd = zeros(n, min(m,n));
for r = 1:min(m,n)
    xtsvd(:,r) = V(:,1:r) * (diag(S(1:r,1:r)).^-1 .* ((U(:,1:r))' * b));
end
