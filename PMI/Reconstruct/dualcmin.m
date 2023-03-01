%DUALCMIN       Dual wavelength regularization: constant ratio
%
%   [x1, x2] = dualcmin(A1, b1, A2, b2, c, lambda)
%
%   x1, x2      The two estimates of the original distribution.
%
%   A1, b1      The system and data for the 1st wavelength.
%
%   A2, b2      The system and data for the 2nd wavelength.
%
%   c           The ratio between x1 and x2.
%
%   lambda      The regularization parameter.
%
%
%   DUALCMIN solves the minimization problem
%                       2                2               2 
%    min   ||A1x1 - b1||  + ||A2x2 - b2||  + ||x1 - cX2||
%   x1,x2

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
%  $Log: dualcmin.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.1  1998/06/03 16:08:40  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, x2] = dualcmin(A1, b1, A2, b2, c, lambda)

[M N] = size(A1);
M1 = A1'*A1 + lambda*eye(N);
M2 = A2'*A2 + c^2*lambda* eye(N);
x1 = (M2*M1 -  c^2 * lambda^2 * eye(N)) \ (M2*A1'*b1 + lambda*c*A2'*b2);

x2 = 1/(lambda*c) * (M1*x1 - A1'*b1);
