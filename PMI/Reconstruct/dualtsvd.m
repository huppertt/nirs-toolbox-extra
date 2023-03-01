%DUALTSVD       2 System truncated SVD solution for a fat matrix problems.
%
%   [x1 x2 rmsdiff] = dualtsvd(A1, A2, b1, b2, c, rmax, 
%                         U1, S1, V1, U2, S2, V2);
%
%   x1          The estimate of the first system.
%
%   x2          The estimate of the second system.
%
%   rmsdiff     The RMS of the difference between the scaled estimates.
%
%   A1          The forward matrix for the first system.
%
%   A2          The forward matrix for the second system.
%
%   b1          The data for the first system.
%
%   b2          The data for the first system.
%
%   c           The ratio x1/x2.
%
%   rmax        The index of the last singular values to use.
%
%   U1,S1,V1    OPTIONAL: The economy SVD's for A1 and A2.
%   U2,S2,V2
%
%   Calls: fattsvd
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
%  $Log: dualtsvd.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1998/04/30 19:09:34  rjg
%  Change difference measure to RMS.
%  
%  Fixed a bug when the SVD was supplied on the command line.
%
%  Revision 1.1  1998/04/30 18:47:53  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1, x2, rmsdiff] = dualtsvd(A1, A2, b1, b2, c, rmax, U1, S1, V1, ...
    U2, S2, V2)

rmsdiff = zeros(rmax, 1);

%%
%%  Compute the economy SVDs of each system if necessary
%%
if nargin < 12
    [x2 U2 S2 V2] = fattsvd(A2, b2, 1);
else
    x2 = fattsvd(A2, b2, 1, U2, S2, V2);
end
if nargin < 9
    [x1 U1 S1 V1] = fattsvd(A1, b1, 1);
else
    x1 = fattsvd(A1, b1, 1, U1, S1, V1);
end
rmsdiff(1) = sqrt(mean((x1 - c * x2).^2));

%%
%%  Calculate all of the truncation values upto r
%%
for iR = 2:rmax,
    x1 = fattsvd(A1, b1, iR, U1, S1, V1);
    x2 = fattsvd(A2, b2, iR, U2, S2, V2);
    rmsdiff(iR) = sqrt(mean((x1 - c * x2).^2));
end
