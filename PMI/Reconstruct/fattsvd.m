%FATTSVD        Truncated SVD solution to a fat matrix problem.
%
%   [xtsvd U S V]= fattsvd(A, b, r, U, S, V)
%
%   xtsvd       The estimate of x (n x 1).
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
%   r           The truncation parameter (1 <= r <= m).
%
%   Calls: none.
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/09/06 12:09:09 $
%
%  $Revision: 1.4 $
%
%  $Log: fattsvd.m,v $
%  Revision 1.4  2000/09/06 12:09:09  dboas
%  Removed fattsvdls.m and changed fattsvd.m to deal with options such as
%  least squares.
%
%  Revision 1.3  2000/07/27 15:11:09  dboas
%  Minor fix to use SVDS properly.
%
%  Revision 1.2  2000/07/21 00:23:25  dboas
%
%  Now using SVDS instead of SVD because of a convergence problem arising from
%  trying to get the entire singular value spectrum.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.3  1999/02/15 19:42:43  rjg
%  Added ability to handle tall and square problems correctly.
%
%  Revision 1.2  1998/12/22 16:55:56  rjg
%  Updated help section.
%
%  Revision 1.1  1998/04/29 17:17:45  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xtsvd, U, S, V, Binv] = fattsvd(A, b, r, option, U, S, V)


%%
%%  Compute the economy SVD of the tall problem if necessary
%%
if nargin < 7
  if option.Lsq
    B = A'*A;
  else
    B = A;
  end
  
  [m n] = size(B);
  if m < n
    if option.FullSVS
      [V S U] = svd(B', 0);
    else
      [V S U] = svds(B', r+1);
    end
  else
    if option.FullSVS
      [U S V] = svd(B, 0);
    else
      [U S V] = svds(B, r+1);
    end
  end
end

Binv = V(:,1:r) * S(1:r,1:r)^-1 * U(:,1:r)';

if option.Lsq
  xtsvd = Binv * (A'*b);
else
  xtsvd = Binv * b;
end
