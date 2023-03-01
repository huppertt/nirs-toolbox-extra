%JTSVD          Joint Truncated SVD solution.
%
%   [xtsvd U S V] = jmtsvd(A, b, r, U, S, V)
%
%   xtsvd       The truncated SVD solutions.
%
%   U, S, V     OPTIONAL: The economy SVD of the block diagonal system.
%
%   A           The two system matrices in a 3D array. A_1 = A(:,:,1),
%               A_2 = A(:,:,2).
%
%   r           The set of truncation parameters to compute.
%
%
%   JTSVD computes the SVD of the block diagonal system unless the SVD
%   is provided.  The estimates are computed from the truncated SVD of
%   the block diagonal system.
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
%  $Log: jtsvd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.0  1998/09/22 17:58:12  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xtsvd, U, S, V] = jtsvd(A, b, r, U, S, V)

m = size(A,1);
n = size(A,2);

%%
%%  Compute the SVD if not supplied.
%%  This only handles underdetermined systems
%%
if nargin < 6
    [V1 S1 U1] = svd(A(:,:,1)', 0);
    [V2 S2 U2] = svd(A(:,:,2)', 0);
    s1 = diag(S1);
    s2 = diag(S2);
    [sr ipr] = sort([s1; s2]);
    iPermute = rev(ipr);

    U = [U1 zeros(size(U1)); zeros(size(U1)) U2];
    U = U(:, iPermute);
    V = [V1 zeros(size(V1)); zeros(size(V1)) V2];
    V = V(:, iPermute);
    S = diag(rev(sr));
end


%%
%%  Compute the truncated SVD solutions for all truncation parameters.
%%
rmax = max(r);
rmin = min(r);
nr = length(r);
disp('Computing truncated SVD solution...');
fc = (U(:,1:rmax)' * b) ./ diag(S(1:rmax,1:rmax));

xtsvd = zeros(2*n, nr);
for i = 1:nr
    xtsvd(:,i) = V(:,1:r(i)) * fc(1:r(i));
end
