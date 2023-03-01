%MNQRSET        Compute the min. norm solutions for a set of matrices.
%
%   xmn = mnqrset(A, b, nColMin, nColMax)
%
%   A           The underdetermined forward matrix.
%
%   b           The family of right hand sides.  The first column is
%               associated with the smallest system (the right most
%               nColMin columns of A).  The last column is the right hand
%               side for the largest system (the right most nColMax
%               columns of A).  If b is a single column it will be
%               employed as the right hand side for all systems.
%
%   nColMin     The minimum number of columns to employ for
%               calculating the minimum norm solution.
%
%   nColMax     OPTIONAL: The maximum number of columns to employ for
%               calculating the minimum norm solution (default: the
%               number of columns in A).
%
%   MNQRSET efficiently computes the minimum norm solution to at set of
%   matrices (and optionally a set of right hand sides).  Each matrix
%   must be created by prepending a column to the previous forward
%   matrix.  The routine computes the first minimum norm solution by
%   computing the QR decomposition of the right most columns of A.  Each
%   succesive QR decompistion is computed adding the new column and using
%   givens rotations to zero the first subdiagonal in R created by
%   prepending the new column.
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
%  $Log: mnqrset.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.1  1999/11/10 18:16:04  rjg
%  Removed use of newR and z temp vars to reduce memory usage.
%
%  Revision 1.0  1998/09/22 17:58:52  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xmn = mnqrset(A, b, nColMin, nColMax)

%%
%%  Initializations
%%
[m n] = size(A);
if nargin < 4
    nColMax = n;
end
r = nColMin:nColMax;
nSoln = length(r);
if size(b,2) ~= nSoln
    b = repmat(b,1, nSoln);
end
xmn = zeros(nColMax, nSoln);

%%
%%  Compute the first min norm solution via the QR decomposition of A'
%%
[Q R] = qr(A(:,n-nColMin+1:n)', 0);
xmn(nColMax-nColMin+1:nColMax, 1) = Q *  (R' \ b(:,1));

%%
%%  Loop over each succesive column of A
%%
for i = 2:nSoln
    %%
    %%  Insert the previous QR decomposition in the newQ and newR
    %%  matrices
    %%
    newQ = zeros(r(i),m);
    newQ(1,1) = 1;
    newQ(2:r(i),2:m+1) = Q;
    Q = newQ;
    clear newQ
    R = [A(:,n-r(i)+1)'; R];

    %%
    %%  Compute the givens rotations to make the new R upper diagonal
    %%
    for j = 1:m
        [c s] = givens(R(j,j), R(j+1,j));
        %%
        %%  Update R
        %%
        r1 = c * R(j,:) - s * R(j+1,:);
        R(j+1,:) = c * R(j+1,:) + s * R(j,:);
        R(j,:) = r1;
        
        %%
        %%  Update Q
        %%
        c1 = c * Q(:,j) - s * Q(:,j+1);
        Q(:,j+1) = c * Q(:,j+1) + s * Q(:,j);
        Q(:,j) = c1;
    end
    
    %%
    %%  Force R to be upper diag for the back slash operator
    %%
    R = triu(R(1:m,:));
    Q = Q(:,1:m);

    %%
    %%  Compute the min. norm solution
    %%
    xmn(nColMax-r(i)+1:nColMax,i) = Q * (R' \ b(:,i));
end


%%
%%  From Golub & Van Loan
%%
function [c,s] = givens(a,b)
if b == 1
    c = 1; s = 0;
else
    if abs(b) > abs(a)
        tau = -a/b;
        s = 1/sqrt(1+tau^2);
        c = s * tau;
    else
        tau = -b/a;
        c = 1/sqrt(1+tau^2);
        s = c * tau;
    end
end
