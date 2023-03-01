%DUALGSVD       Compute the GSVD solution to a regularized, dual wavelength.
%
%   [x U V X C S] = dualgsvd(A, L, b, lambda, X, C, S)
%
%   x           The estimate of x.
%
%   U,V,X,C,S   The gsvd of A,L.
%
%   lambda      The value(s) of the regularization parameters for which to
%               solve the minimization problem.
%
%
%   DUALGSVD attempts to solve the minimization problem.
%
%   min ||Ax - b|| + \lambda * ||Lx||
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
%  $Log: dualgsvd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.0  1998/09/22 17:54:48  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, U, V, Xinv, C, S] = dualgsvd(A, L, b, lambda, Xinv, C, S)

[m p] = size(A);
[n pl] = size(L);
if pl ~= p
    error('A and L must have the same number of columns');
end

if nargin < 7
    [U V X C S] = gsvd(A, L, 0);
    Xinv = inv(X);
end

nLambda = length(lambda);
x = zeros(p, nLambda);
Rterm = Xinv * (A' * b);

%%
%%  This needs to be more robust than just diag, also what about when there
%%  is no inverse to X.
%%
c_vec = (diag(C)).^2;
s_vec = (diag(S)).^2;
if length(s_vec) < p
    s_vec(p) = 0;
end
lamba = lambda .^2;

for iLambda = 1:nLambda,
    x(:,iLambda) = Xinv' * ((c_vec + lambda(iLambda)*s_vec).^-1 .* Rterm);
end

