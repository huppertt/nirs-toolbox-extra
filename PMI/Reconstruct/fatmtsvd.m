%FATMTSVD       Modified truncated SVD solution with null space wieghting.
%
%   [xlm, xtsvd, xns, U, S, V] = fatmtsvd(A, b, r, L, lambda, U, S, V)
%
%   xlm         The minimized semi-norm solution.
%
%   xtsvd       The truncated SVD solution.
%
%   xns         The null space vector subtracted from the truncated SVD
%               solution meet the constriant.
%
%   U,S,V       The SVD of A.
%
%   A,b         The system to be solved.
%
%   r           The number of singular values/singular vectos to use for the
%               truncated SVD system.
%
%   L           The semi-norm operator.
%
%   lambda      A parameter that controls how much of the null space is
%               added to the truncated SVD solution.
%
%   Reference: Piecewise Polynomial Solutions Without a priori Break Points
%              P.C. Hansen, Klaus Mosegaard,  Numerical Linear Algebra with
%              Applications, Vol 3(6) 513-524 1996
%
%
%   Calls: none.
%
%   Bugs: not sure if this solves the short-fat problem correctly?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: fatmtsvd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1998/09/22 17:52:30  rjg
%  Now handles cases where the null space coefficient calculation requires a
%  least squares approach as opposed to a min norm approach.
%
%  Revision 1.1  1998/06/03 16:12:22  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xlm, xtsvd, xns, U, S, V] = fatmtsvd(A, b, r, L, lambda, U, S, V)

[M N] = size(A);
if nargin < 5
    lambda = 1;
end

%%
%%  Compute the SVD of the forward matrix
%%
if nargin < 8
    disp('Computing full SVD...');
    [U S V] = svd(A);
end

%%
%%  Compute the minimum norm solution to the truncated system.
%%
disp('Computing truncated SVD solution...');
xtsvd = V(:,1:r)* (diag(S(1:r,1:r)).^-1 .* ((U(:,1:r))' * b));

%%
%%  Check to see if the semi-nomrm problem is underdetermined.
%%  If so compute the minimim norm solution for y, otherwise compute
%%  the least squares solution.
%%
Lm = size(L,1);
if Lm < N-r
    disp('Computing the minimum norm solution for the null space...');
    y = fatmn(L*V(:,r+1:N), L*xtsvd);
else
    disp('Computing the least squares solution for the null space...');
    y = (L * V(:,r+1:N)) \ (L * xtsvd);
end

%%
%%  Compute the null space component
%%
disp('Adding in the null space component..');
xns = V(:,r+1:N) * y;
xlm = xtsvd - lambda * xns;
