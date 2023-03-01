%ART            Algebraic Reconstruction Technique
%
%   x = art(A, b, xo, nIter, w, iRow)
%
%   x           The estimate of the x vector.
%
%   A           The forward matrix.
%
%   b           The measured data.
%
%   xinit       OPTIONAL: An initial guess, if not supplied then 0 will be
%               used.
%
%   nIter       OPTIONAL: The maximum number of iterations to compute
%               (default 10 * number of rows).  If nIter is a vector and
%               estimate is returned for each element in nIter.  The number
%               of iterations must be increasing.
%
%   w           OPTIONAL: The relaxation parameter, defining how far to "step"
%               on each iteration.  A value of 1 causes each step to reach
%               the hyperplane of the current orthogonal projection.  Less than
%               1 cause the step to fall short of the hyperplane.  The default
%               value is 1.
%
%   iRow        OPTIONAL: The initial row to project onto.  Useful for
%               examining convergence performance.
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
%  $Log: art.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.5  1999/11/10 18:19:44  rjg
%  Added constant relaxation parameter w.
%
%  Revision 1.4  1999/11/05 20:24:56  rjg
%  Transpose A so that the rows of A are in cache (sequential in memory)
%  instead of the columns.
%
%  Revision 1.3  1999/01/04 22:00:35  rjg
%  Fixed help text.
%
%  Revision 1.2  1998/09/14 18:43:03  rjg
%  Added the ability to process many iteration estimates.
%
%  Revision 1.1  1998/06/03 16:08:15  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = art(A, b, xinit, nIter, w, CurrRow)

[nr nc] = size(A);
nEst = length(nIter);

%%
%%  Default arguments
%%
if nargin < 6
    CurrRow = 1;
    if nargin < 5
        w = 1;
        if nargin < 4
            nIter = 10 * nr;
            if nargin < 3
                xinit = zeros(nc, 1);
            end
        end
    end
end
x = zeros(nc, nEst);

%%
%%  Transpose A so that all of the operations are now column operations, should
%%  work in cache more
%%
A = A.';

%%
%%  Precompute the row norm
%%
rownorm = zeros(nr,1);
for i=1:nr,
    rownorm(i) = A(:,i).' * A(:,i);
end

%%
%%  Loop over the number of iterations requested
%%
for j = 1:nEst
    if j == 1
        x(:, j) = xinit;
        N = nIter(1);
    else
        x(:,j) = x(:,j-1);
        N = nIter(j) - nIter(j-1);
    end

    %%
    %%  Iterate over each row 
    %%
    for i = 1:N
        x(:,j) = x(:,j) + w * (b(CurrRow) - x(:,j).' * A(:, CurrRow)) / ...
            rownorm(CurrRow) * A(:, CurrRow);
        CurrRow = CurrRow + 1;
        if CurrRow > nr
            CurrRow = 1;
        end
    end
end
