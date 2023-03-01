%SIRT           Simultaneous iterative reconstruction technique
%
%   x = sirt(A, b, xo, nIter, w)
%
%   x           The estimate of the x vector.
%
%   A           The forward matrix.
%
%   b           The measured data.
%
%   xo          OPTIONAL: An initial guess, if not supplied then 0 will be
%               used.
%
%   nIter       OPTIONAL: The maximum number of iterations to compute
%               (default 10 * number of rows).
%
%   w           OPTIONAL: The relaxation parameter, defining how far to "step"
%               on each iteration.  A value of 1 causes each step to reach
%               the hyperplane of the current orthogonal projection.  Less than
%               1 cause the step to fall short of the hyperplane.  The default
%               value is 1.
%
%
%   SIRT
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
%  $Log: sirt.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.6  1999/11/18 14:19:56  tgaudett
%  Deal with zero row norm
%
%
%	Revision 1.6  1999/10/12	DAB
%	Deal with zero signal rows
%
%  Revision 1.5  1999/11/11 00:10:58  rjg
%  Moved x(:,j) out of the inside for loop.
%
%  Revision 1.4  1999/11/10 20:18:14  rjg
%  Added relaxation parameter.
%
%  Revision 1.3  1998/09/14 18:43:38  rjg
%  Added the ability to compute many iteration estimates.
%
%  Revision 1.2  1998/06/03 16:15:33  rjg
%  Changed while loop to for loop.
%  Combined current x assignement and division by the number of rows for
%  the average step direction inside the loop.
%
%  Revision 1.1  1998/04/29 15:56:43  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = sirt(A, b, xinit, nIter, w)

[nr nc] = size(A);
nEst = length(nIter);

%%
%%  Default arguments
%%
if nargin < 5
    w = 1;
    if nargin < 4
        nIter = 10 * nr;
        if nargin < 3
            xinit = zero(nc, 1);
        end
    end
end
x = zeros(nc, nEst);

%%
%%  Transpose A so that all of the operations are now column operations, should
%%  work in cache more
A = A.';

%%
%%  Precompute the row norm
%%
rownorm = zeros(nr,1);
for i=1:nr,
   rownorm(i) = A(:,i).' * A(:,i);
  	if rownorm(i)==0 
		rownorm(i) = 9e50;	%% deal with zero signal rows - DAB 99-10-12
	end

end

%%
%%  Loop over the number of iterations requested
%%
for j = 1:nEst
    %%
    %%  Copy the estimate from the previous iteration parameter to the ...
    %%  estimate for the new x(:,j)iteration parameter.
    %%
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
        RelResid = (b - (x(:,j).' * A).') ./ rownorm;
        xnew = zeros(nc, 1);
        for iRow = 1:nr
            xnew = xnew +  RelResid(iRow) * A(:,iRow);
        end
        x(:,j) = x(:,j) + w / nr * xnew;
    end
end
