%TSVDCEST       Estimate the constant ratio between two vectors by the TSVD.
%
%   [c_est c_std U S V] = tsvdcest(A, b, k, thresh, U, S, V)
%
%   c_est       The estimated values of C.
%
%   c_std       The standard deviation of each of the estimates of C.
%
%   U S V       OPTIONAL: The economy SVD of each of the forward matrices.
%
%   A           A m x n x 2 array containing the two forward matrices.
%
%   b           A m x 2 array containing the two measurement vectors.
%
%   k           A vector of truncation parameters to compute.
%
%   thresh      OPTIONAL: The threshold of the peak reconstruction value
%               for detection of an object (default: 0.5).
%
%
%   TSVDCEST
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
%  $Log: tsvdcest.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.0  1998/09/22 18:04:21  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c_est, c_std, U, S, V] = tsvdcest(A, b, k, thresh, U, S, V)

nK = length(k);
c_est = zeros(nK,1);
c_std = zeros(nK,1);

%%
%%  Compute each of the SVDs if necessary.
%%
if nargin < 5
    for i = 1:2
        [V(:,:,i) S(:,:,i) U(:,:,i)] = svd(A(:,:,i)', 0);
    end
end

%%
%%  Loop over each of the truncation parameters
%%
for i = 1:length(k)

    %%
    %%  Compute the tsvd soln for each wavelength
    %%
    xtsvd1 = V(:,1:k(i),1) * (diag(S(1:k(i),1:k(i),1)).^-1 .* ...
        ((U(:,1:k(i),1))' * b(:,1)));
    xtsvd2 = V(:,1:k(i),2) * (diag(S(1:k(i),1:k(i),2)).^-1 .* ...
        ((U(:,1:k(i),2))' * b(:,2)));
        
    %%
    %%  Threshold both reconstruction vector
    %%
    peak1 = max(xtsvd1);
    peak2 = max(xtsvd2);
    idx = (xtsvd1 > thresh*peak1) & (xtsvd2 > thresh*peak2);
    %%
    %%  Compute the mean ratio between the two thresholded regions
    %%  as well as the the standard deviation
    %%
    if sum(idx) > 0
        ratio = xtsvd1(idx) ./ xtsvd2(idx);
        c_est(i) = mean(ratio);
        c_std(i) = std(ratio);
    end
end
