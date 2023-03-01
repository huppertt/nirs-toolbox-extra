%FATMN          Minimum norm solution to a fat system using a direct solution.
%
%   xmn = fatmn(A, b)
%
%   xmn         The minimum norm solution.
%
%   A           The under-determined matrix (short and fat, m > n).
%
%   b           The measurement vector.
%
%
%   FATMN computes the minimum norm solution to the underdetermined system Ax=b
%   by solving the equation
%
%       xmn = A' * (A * A')^-1 * b
%
%   The inverse is not actually computed.  The associated problem
%
%       (A * A') * z = b    =>  z = (A * A')^-1 * b
%
%   is solved for z by Gaussian elimination.  The min. norm solution is then
%   given by
%
%       xmn = A' * z
%
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
%  $Log: fatmn.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.1  1998/04/29 17:17:10  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xmn = fatmnqr(A, b)
z = (A * A') \ b;
xmn = A' * z;

