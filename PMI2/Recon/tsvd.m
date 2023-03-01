% TSVD   Truncated singular value decomposition
%
% [X,U,S,V] = tsvd(A, Y, nSV, [u, s, v]);
%
% Y     - residues
% A     - forward matrix (not used if u,s,v are passed)
% nSV   - number of singular values to use (can be a vector)
% u,s,v - OPTIONAL: cached SVD of A
%
% X     - reconstructed image
% U,S,V - sigular values used in reconstruction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang,
%                     Jonathan Stott
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[X,U,S,V] = tsvd(A, Y, nSV, u, s, v);

% Sanity check initial arguments

if (size(Y,1) == 1 & size(Y,2) ~= 1)
   % ' or .' ???
   Y = Y.';
end

% Fill in remaining arguments as needed

if (exist('u','var') & exist('s','var') & exist('v','var'))
   haveUSV = 1;
else
   haveUSV = 0;
   
   if (~exist('u','var'))
      u = [];
   end
   
   if (~exist('s','var'))
      s = [];
   end
   
   if (~exist('v','var'))
      v = [];
   end
   
   if (~isempty(u) | ~isempty(s) | ~isempty(v))
      error('Either none of or all three of u,s,v must be defined');
   end
end

% Invert the data using truncated SVD

mxSV = max(nSV);
   
if (haveUSV)
   if (mxSV > size(u,2))
      error('Asked for more eigenvectors than passed to function');
   else
      u = u(:,1:mxSV);
      s = s(1:mxSV,1:mxSV);
      v = v(:,1:mxSV);
   end
else
   if (mxSV <= 20)
      [u,s,v] = svds(A,mxSV);
   else
      [u,s,v] = svd(A);
      
      u = u(:,1:mxSV);
      s = s(1:mxSV,1:mxSV);
      v = v(:,1:mxSV);
   end
end

for n = 1:length(nSV)
   % t is sparse, so having extra elements won't slow down the multiplications
   t = 0*speye(mxSV);
   
   for k = 1:nSV(k)
      t(k,k) = 1 ./ s(k,k);
   end

   X(:,n) = v * t' * u' * Y;
end

U = u;
S = sparse(s);
V = v;

return;



