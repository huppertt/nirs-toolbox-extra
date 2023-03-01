%DUALREG        Dual lambda regularized solution
%
%   [x1, x2] = dualreg(A1, b1, A2, b2, c, l1, l2, L)
%
%   x1, x2      The two estimates of the original distribution.
%
%   A1, b1      The system and data for the 1st wavelength.
%
%   A2, b2      The system and data for the 2nd wavelength.
%
%   c           The ratio between x1 and x2.
%
%   l1          The regularization parameter controlling the scaled
%               difference between the two solutions.
%
%   l2          The regularization parameter controlling the seminorm of x1.
%
%   L           The semi-norm operator.
%
%   DUALREG solves the minimization problem
%                       2                2                   2             2
%    min   ||A1x1 - b1||  + ||A2x2 - b2||  + l_1 ||x1 - cx2|| + l_2 ||Lx1||
%   x1,x2
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
%  $Log: dualreg.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.0  1998/09/22 17:55:45  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, x2] = dualreg(A1, b1, A2, b2, c, l1, l2, L)

%%
%%  Precompute the normal equations for each system.
%%
A1n = A1' * A1;
A2n = A2' * A2;
Ln = L' * L;

%%
%%  Compute the system that specifies x1
%%
LHS = A2n * A1n + l2 * A2n * Ln + l1 * A2n + (l1*c^2) * A1n + (l1*l2*c^2) * ...
    Ln;

RHS = l1*c*(A2' * b2) + A2n * (A1' * b1) + l1 * c^2 * (A1' * b1);
x1 = LHS \ RHS;

%%
%%  Compute x2 from x1
%%
x2 = (A1' * b1 - A1n * x1 + l2 * Ln * x1 + l1 * x1) ./ (l1 * c);
