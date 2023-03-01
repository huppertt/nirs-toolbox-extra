%JCTSVD         Joint constrained TSVD solution.
%
%   [xctsvd U S V] = jctsvd(A, b, r, c, U, S, V)
%
%   xctsvd       The truncated constrianed SVD solutions.
%
%   U, S, V     OPTIONAL: The economy SVD of the block diagonal system.
%
%   A           The two system matrices in a 3D array. A_1 = A(:,:,1),
%               A_2 = A(:,:,2).
%
%   r           The set of truncation parameters to compute.
%
%   c           The expected ratio between x_1 and x_2 as x_1 / x_2.
%
%
%   JCTSVD computes the TSVD solution to
%
%   |A1  0| |x1|   |b1|
%   |0  A2| |x2| = |b2|
%   |I -cI|        |0 |
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
%  $Log: jctsvd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.0  1998/09/22 17:56:39  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [xctsvd, U, S, V] = jmtsvd(A, b, r, c, U, S, V)

m = size(A, 1);
n = size(A, 2);
b = [b(:)
     zeros(n,1) ];

%%
%%  Compute the SVD of the large system if necessary
%%
if nargin < 7
    ts_svd = clock;
    A = [diagdsa(A)
        eye(n) -c*eye(n)];
    %%
    %%  Compute the most efficient economy SVD
    %%
    if 2*n > n+2*m
        disp('Underdetermined')
        [V S U] = svd(A', 0);
    else
        disp('Overdetermined')
        [U S V] = svd(A, 0);
    end
    te_svd = clock;
end

%%
%%  Compute the truncated SVD solutions for all truncation parameters.
%%
rmax = max(r);
rmin = min(r);
nr = length(r);
disp('Computing truncated SVD solution...');
fc = (U(:,1:rmax)' * b) ./ diag(S(1:rmax,1:rmax));

xctsvd = zeros(2*n, nr);
for i = 1:nr
    xctsvd(:,i) = V(:,1:r(i)) * fc(1:r(i));
end
te_soln = clock;
fprintf('SVD time %f\n', etime(te_svd, ts_svd));
fprintf('SOLN time %f\n', etime(te_soln, te_svd));
