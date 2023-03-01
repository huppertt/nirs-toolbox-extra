% SHOWIMAGE  Show reconstructed data sets
%
% h = showimage(Medium, [Mua], [Musp]);
%
% Medium   - PMI Structure
% Mua,Musp - reconstructed data.  One of the two must be passed
% h        - figure handle of new figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2004, David Boas, Dana Brooks, Rick Gaudette, 
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

function[fig] = showImage(Medium, Mua, Musp);

[X, Y, Z] = sampleVolume(Medium.CompVol);
uX = unique(X);
uY = unique(Y);
uZ = unique(Z);

% System size

nx = length(uX);
ny = length(uY);
nz = length(uZ);
nv = length(X);

clear X Y Z ;

% What do we need to plot

doMua = 0;
doMus = 0;

if (exist('Mua','var') & (~isempty(Mua)))
   doMua = 1;
   Mua   = Mua(:);

   if (length(Mua) ~= nv & length(Mua) ~= 2*nv)
     error('showImage() does not support irregular grid spacings');
   end
end

if (exist('Musp','var') & (~isempty(Musp)))
   doMus = 1;
   Musp  = Musp(:);

   if (length(Musp) ~= nv & length(Musp) ~= 2*nv)
     error('showImage() does not support irregular grid spacings');
   end
end

if (doMua & ~doMus & length(Mua) == 2*nv)
   % Mus and Mua have been packed together into a single vector
   Musp = Mua(nv + [1:nv]);
   Mua  = Mua([1:nv]);
   
   doMus = 1;
end

if (doMus & ~doMua & length(Musp) == 2*nv)
   % Mus and Mua have been packed together into a single vector
   Mua  = Musp([1:nv]);
   Musp = Musp(nv + [1:nv]);
   
   doMua = 1;
end

if (~doMua & ~doMus)
   error('Empty data sets; I have no work to do!');
end

% Size of the subplots

[M,N] = plotDimensions(nz, doMua, doMus);

% Open the figure;

fig = figure;

% Plot the data

if (doMua)
   Mua = reshape(Mua, ny, nx, nz);
   
   plotSet(uX, uY, uZ, M, N, 1:nz, Mua, 1);
end

if (doMus)
   if (doMua)
      offset = N*(M/2);
   else
      offset = 0;
   end
   
   Musp = reshape(Musp, ny, nx, nz);

   plotSet(uX, uY, uZ, M, N, offset + [1:nz], Musp, 0);
end

drawnow;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSet(uX, uY, uZ, M, N, zlst, Mu, isAbs)

maxMu = max(Mu(:));
minMu = min(Mu(:));

if (minMu == maxMu)
   % Try to plot constants a little more intellegently
   minMu = 0.9 * minMu;
   maxMu = 1.1 * maxMu;
end

for k = 1:length(zlst)
   z = zlst(k);
   zp = z - zlst(1) + 1;
   
   subplot(M,N,z);
      
   imagesc(uX, uY, Mu(:,:,zp), [minMu maxMu]);
   set(gca,'ydir','normal');
   
   if (isAbs)
      title(sprintf('\\mu_a  Z=%5.3f', uZ(zp)));
   else
      title(sprintf('\\mu_s'' Z=%5.3f', uZ(zp)));
   end
   
%  colorbar;
end
   
colorbar;

return;
