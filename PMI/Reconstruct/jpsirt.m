%JPSIRT         Joint Parallel SIRT solution for two wavelengths.
%
%   [x c]= jsirt(A, b, xo, nIter, nBlk, C)
%
%   x           The estimate of the x vector.
%
%   C
%
%   A,b         The joint underdetermined system of equations to be solved.
%               A has the form A_1 = A(:,:,1), A_2 = A(:,:,2). b should contain
%               
%
%   xo          OPTIONAL: An initial guess, if not supplied then 0 will be
%               used.
%
%   nIter       OPTIONAL: The number of interations to compute for each system
%               before switching to the next system (default 10 * # of rows).
%
%   n
%
%   JPSIRT
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
%  $Log: jpsirt.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.0  1998/09/22 17:57:17  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, C, xmean] = jpsirt(A, b, x, nIter, nBlk, C)


m = size(A,1);
n = size(A,2);

%%
%%  Precompute the row norm
%%
rownorm1 = zeros(m,1);
for i=1:m
    rownorm1(i) = A(i,:,1) * A(i,:,1).';
end
rownorm2 = zeros(m,1);
for i=1:m
    rownorm2(i) = A(i,:,2) * A(i,:,2).';
end


x1 = x(:,1);
x2 = x(:,2);
for i = 1:nBlk
    %%
    %%  SIRT on the two systems
    %%
    x1 = sirt_iter(A(:,:,1), b(:,1), x1, rownorm1, nIter);
    x2 = sirt_iter(A(:,:,2), b(:,2), x2, rownorm2, nIter);

    %%
    %%  Estimate of C, and then average the two results
    %%
    if nargin < 6
        [x1 x2 xmean C] = cest(x1, x2);
    else
        [x1 x2 xmean] = cest(x1, x2, C);
    end
end

x = [x1 x2];
return

%%
%%   SIRT function
%%
function x = sirt_iter(A, b, x, rownorm, nIter)
[m n] = size(A);
for j = 1:nIter
    RelResid = (A*x - b) ./ rownorm;
    xnew = zeros(n, 1);
    for iRow = 1:m
        xnew = xnew + (x - RelResid(iRow) * A(iRow,:).');
    end
    x = xnew ./ m;
end

return

%%
%%  Constant ratio estimate and averaging of the two estimates
%%
function [x1, x2, C, xmean] = cest(x1, x2, C)

%%
%%  Estimate of constant ratio
%%
if nargin < 3
    peak1 = max(x1);
    peak2 = max(x2);

    idx = (x1 > 0.5*peak1) & (x2 > 0.5*peak2);
    C = mean(x1(idx) ./ x2(idx))
end

%%
%%  Compute the estimate of the two systems from their average
%%
xmean = mean([x1'; x2'])';
x1 = 2 * xmean ./ (1 + 1/C);
x2 = (1/C) * x1;
return
