%JMTSVD         Joint MTSVD estimate of a pair of systems.
%
%   [xlm, xtsvd, U, S, V] = jmtsvd(A, b, r, c, U, S, V)
%
%   xlm         The Joint TMSVD estimate(s).
%
%   xtsvd       The truncated SVD estimate(s).
%
%   U S V       OPTIONAL: The SVD of the diagonal system.
%
%   A,b         The joint underdetermined system of equations to be solved.
%               A has the form A_1 = A(:,:,1), A_2 = A(:,:,2).  b should either
%               be a stacked vector of the two RHS or a 2 column matrix.
%
%   r           The truncation parameters to compute.
%
%   c           The expected ratio between x_1 and x_2 as x_1 / x_2.
%
%   JMTSVD computes the Joint MTSVD estimate of pair of systems of equations
%   with object functions related by a constant value (c).
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
%  $Log: jmtsvd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1999/11/10 17:47:21  rjg
%  Added intermediate clean up of temporary variables to save space.
%
%  Revision 1.1  1998/09/22 17:45:10  rjg
%  Computing of the null space components for multiple truncation values is
%  now much more efficent.  mnqrset computes a single full QR decomposition
%  and appropriately updates the QR as more of the null space is used.
%
%  Revision 1.0  1998/09/17 19:46:02  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xlm, xtsvd, U, S, V] = jmtsvd(A, b, r, c, U, S, V)

b = b(:);
m = size(A,1);
n = size(A,2);
rmax = max(r);
rmin = min(r);
nr = length(r);

%%
%%  Compute the SVD if not supplied.
%%  This only handles underdetermined systems
%%
ts_svd = clock;
if nargin < 7
    disp('Computing SVDs...');
    nNS = n - m;
    [U1 S1 V1] = svd(A(:,:,1));
    [U2 S2 V2] = svd(A(:,:,2));
    s1 = diag(S1);
    s2 = diag(S2);
    [sr ipr] = sort([s1; zeros(nNS,1); s2; zeros(nNS,1)]);
    iPermuteR = rev(ipr);
    [sr ipr] = sort([s1; s2]);
    iPermuteL = rev(ipr);

    U = [U1 zeros(size(U1)); zeros(size(U1)) U2];
    U = U(:, iPermuteL);
    V = [V1 zeros(size(V1)); zeros(size(V1)) V2];
    V = V(:, iPermuteR);
    S = [diag(rev(sr)) zeros(2*m,2*(n-m))];
end


%%
%%  Compute the truncated SVD solutions for all truncation parameters.
%%
disp('Computing truncated SVD solution...');
fc = (U(:,1:rmax)' * b) ./ diag(S(1:rmax,1:rmax));

xtsvd = zeros(2*n, nr); 
for i = 1:nr
    xtsvd(:,i) = V(:,1:r(i)) * fc(1:r(i));
end
clear fc iPermuteR iPermuteL sr ipr s2 s1 U2 S2 V2 U1 S1 V1 A
te_svd = clock;
if nargout < 3
    clear U S
end

%%
%%  Compute the minimum norm solution to the weighted system for
%%  calculating the null space component
%%
disp('Computing null space weights...');
xmn = mnqrset(V(1:n,rmin+1:2*n) - c*V(n+1:2*n,rmin+1:2*n), ...
    fliplr(xtsvd(1:n,:) - c*xtsvd(n+1:2*n,:)), 2*n-rmax, 2*n-rmin);
xmn = fliplr(xmn);
xlm = zeros(2*n, nr);
for i = 1:nr
    xns = V(:,r(i)+1:2*n) * xmn(r(i)-rmin+1:2*n-rmin, i);
    xlm(:,i) = xtsvd(:, i) - xns;
end
te_ns = clock;

fprintf('TSVD time %f\n', etime(te_svd, ts_svd));
fprintf('null space time %f\n', etime(te_ns, te_svd));
