%POSTCONSTRAIN  Post constrain two Slab Image measurement sets.
%
%   xpc = postcontrain(x1, x2, c)
%
%   xpc         The post constrained results.
%
%   x1, x2      The input results to generate the constraint.
%
%   c           The constant ratio to enforce between the sets of results.
%
%
%   x1 = c x2
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
%  $Log: postconstrain.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1999/11/10 18:17:30  rjg
%  Removed DOS linefeeds
%
%  Revision 1.1  1999/09/29 18:39:19  rjg
%  Fixed bug in mean calculation.  Parenthesis where not around x1 and x2!
%
%  Revision 1.0  1998/09/22 17:59:31  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xpc1, xpc2] = postconstrain(x1, x2, c)
[m1 n1] = size(x1);
[m2 n2] = size(x2);

if m1 ~= m2
    error('x1 and x2 are estimates of the same size vector')
end
if n1 ~= n2
    error('x1 and x2 must have the same number of estimates (columns)')
end

xmean = 0.5 * (x1 + x2);
xpc1 = 2/(1+1/c) * xmean;
xpc2 = 1/c * xpc1;
