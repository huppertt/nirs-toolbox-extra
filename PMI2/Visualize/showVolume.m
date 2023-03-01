% SHOWVOLUME  Interactively display orthogonal slices of a 3D volume
%
% h = showVolume(Medium, data);
%
% Medium is the PMI data structure
% data   is a set of datapoints whos size is size specified 
%        by Medium.CompVol.  Typically, this is the reconstructed data.
% h      is the figure handle.
%
% Click in the lower-right plot to end the interactive display

% file hsi3.m by Chuck DiMarzio and Manda Kashambala
%                Northeastern University, December 2000
% file hsi5.m by Jonathan Stott, August 2002
% Renamed showVolume.m and integrated into PMI toolbox, 2003

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

function[h] = showVolume(Medium, data);

[l1 l2 l3 dV] = sampleVolume(Medium.CompVol);
l1 = unique(l1);
l2 = unique(l2);
l3 = unique(l3);

n1 = length(l1);
n2 = length(l2);
n3 = length(l3);

if (length(data(:)) ~= n1*n2*n3)
  error('Size of data vector does not match Medium.CompVol');
end

if (~isreal(data))
  % What about int16 and the like?
  error('Complex data not supported');
end

% data is stored in memory Y-X-Z
data = reshape(data, n2, n1, n3);

lolim = double(min(data(:)));
hilim = double(max(data(:)));

x1 = ceil(n1/2);
x2 = ceil(n2/2);
x3 = ceil(n3/2);

%% Main loop

h = figure;
set(h,'DoubleBuffer','On');

while (1)
   %% Replot the data using the given [x1, x2, x3]
 
   [himg, hcol, hrow, hpix] = plotdata(data, lolim, hilim, ...
				       l1, l2, l3, x1, x2, x3);

   %% Get the new data point
   
   [hin,vin] = ginput(1);
   
   h1 = gca;

   switch (h1)
      case himg
	 x1p = findindex(l1, hin);
	 x2p = findindex(l2, vin);
	 x3p = x3;
      case hcol
	 x1p = x1;
	 x2p = findindex(l2, vin);
	 x3p = findindex(l3, hin);
      case hrow
	 x1p = findindex(l1, hin);
	 x2p = x2;
	 x3p = findindex(l3, vin);
      case hpix
	 % Mouse-click over the spectrum breaks out of the while() loop
	 break;
      otherwise
	 error('Not an image axis');
   end
   
   if     ((x1p < 1) | (x1p > n1))
      disp('X coordinate out of range, using old values');
      disp(x1p);
   elseif ((x2p < 1) | (x2p > n2))
      disp('Y coordinate out of range, using old values');
      disp(x2p);
   elseif ((x3p < 1) | (x3p > n3))
      disp('Z coordinate out of range, using old values');
      disp(x3p);
   else
      x1 = x1p;
      x2 = x2p;
      x3 = x3p;
   end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given a mouse click, find the entry in the list closest to the 
%  given coordinate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[nidx] = findindex(lst, X)

v1 = find(lst < X);
v2 = find(lst > X);

if (isempty(v1) & isempty(v2))
   error('Impossible condition');
end

if      (isempty(v1))
   d1 = lst(2)-lst(1);
   
   if (X < lst(1) - d1)
      % Off scale
      nidx = 0;
   else
      % First voxel still
      nidx = 1;
   end
elseif (isempty(v2))
   d2 = lst(end) - lst(end-1);
      
   if (X > lst(end) + d2)
      % Off scale
      nidx = length(lst) + 1;
   else
      % Last voxel still
      nidx = length(lst);
   end
else
   n1 = v1(end);
   n2 = v2(1);

   % Take whichever index (n1 or n2) is closer to the actual value

   if (X > (lst(n1) + lst(n2))/2)
      nidx = n2;
   else
      nidx = n1;
   end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the set of four images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [himg, hcol, hrow, hpix] = plotdata(data,lolim,hilim,...
					     l1,l2,l3,xi1,xi2,xi3)

x1 = l1(1);   x2 = l1(end);
y1 = l2(1);   y2 = l2(end);
z1 = l3(1);   z2 = l3(end);

himg = subplot(2,2,1);
imagesc(l1, l2, data(:,:,xi3), [lolim hilim]);
axis image
set(gca,'YDir','Normal');
xlabel('X'); ylabel('Y');
title(sprintf('HSI5: Z = %d', xi3));
crosshair(x1, l1(xi1), x2, y1, l2(xi2), y2);
% colorbar;

hcol=subplot(2,2,2);
imagesc(l3, l2, squeeze(data(:,xi1,:)), [lolim hilim]);
axis image;
set(gca,'YDir','Normal');
set(gca,'YAxisLocation','Right');
xlabel('Z'); ylabel('Y');
title(sprintf('X = %d', xi1));
crosshair(z1, l3(xi3), z2, y1, l2(xi2), y2);
   
hrow=subplot(2,2,3);
imagesc(l1, l3, squeeze(data(xi2,:,:))', [lolim hilim]);
axis image
set(gca,'YDir','Normal');
xlabel('X'); ylabel('Z');
title(sprintf('Y = %d',xi2));
crosshair(x1, l1(xi1), x2, z1, l3(xi3), z2);

hpix=subplot(2,2,4);
plot(l3, squeeze(data(xi2,xi1,:)));
axis([l3(1) l3(end) lolim hilim]);
xlabel('Z');
title(sprintf('[X,Y] = [%d,%d]', xi1, xi2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function crosshair(u1,ui,u2,v1,vi,v2)

hold on;

plot([u1,u2],[vi,vi], 'w');
plot([ui,ui],[v1,v2], 'w');

hold off;



