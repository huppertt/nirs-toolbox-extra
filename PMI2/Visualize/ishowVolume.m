% ISHOWVOLUME  Interactively display orthogonal slices of a 3D volume
%
% h = ishowVolume(Medium, data);
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

function[h] = ishowVolume(Medium, data);

% Get parameters necessary to generate the plot

[l1 l2 l3 dV] = sampleVolume(Medium.CompVol);
l1 = unique(l1);
l2 = unique(l2);
l3 = unique(l3);

n1 = length(l1);
n2 = length(l2);
n3 = length(l3);

% Check sizes for self-consistancy

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

% Initial coordinates are in the center of the image

x1 = floor(n1/2) + 1;
x2 = floor(n2/2) + 1;
x3 = floor(n3/2) + 1;

% Set up the figure

h = figure('DoubleBuffer',  'On',  ...
	   'Interruptible', 'Off', ...
	   'BusyAction',    'queue');

udata = cell(1,1);

udata{1} = data;
udata{2} = [ lolim hilim ];
udata{3} = l1;
udata{4} = l2;
udata{5} = l3;
udata{6} = [ x1 x2 x3 ];
udata{7} = [ ];			% Fill in below
udata{8} = @plotdata;
udata{9} = 0;                   % Mid-drag, used by callbacks
udata{10}= @showVolumeCallback;

% Set up the basic plot.  Also assigns the subplot axes to udata{7}.

udata = plotdata(udata);

% Copy back to the figure, DO THIS BEFORE SETTING CALLBACKS!

set(gcf, 'UserData', udata);

% Set callbacks.  Callbacks will operate until window is closed

set(gcf, 'WindowButtonDownFcn', ...
	 'ud = get(gcf,''UserData''); feval(ud{10},1)', ...
         'WindowButtonMotionFcn', ...
	 'ud = get(gcf,''UserData''); feval(ud{10},2)', ...
         'WindowButtonUpFcn',     ...
	 'ud = get(gcf,''UserData''); feval(ud{10},3)');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generic callback routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function showVolumeCallback(mode)

udata = get(gcf, 'UserData');

if (isempty(udata))
   disp('empty udata - returning');
   return;
end

switch (mode)
   case 1
      udata{9} = 1;
      set(gcf, 'UserData', udata, 'Pointer', 'crosshair');
   case 2
      if (udata{9} ~= 0)
         % Button has been pressed, update screen with new [x,y] values

	 udata = updateScreen(udata);
	 set(gcf, 'UserData', udata);
      end
   case 3
      if (udata{9} ~= 0)
         % Record final value (make button click w/o motion work)
         udata = updateScreen(udata);
      end

      udata{9} = 0;
      set(gcf, 'UserData', udata, 'Pointer', 'arrow');
   case 4
      % Close down figure (delete user data);
      set(gcf, 'WindowButtonDownFcn',   '', ...
               'WindowButtonMotionFcn', '', ...
               'WindowButtonUpFcn',     '', ...
               'CloseRequestFcn',       '');

      set(gcf, 'UserData', []);

   otherwise
      error('Unknown callback mode');
      return;
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the set of four images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[udata] = plotdata(Udata)

if (~iscell(Udata))
   error('plotdata passed non-cell argument');
else
   % Copy most parameters across
   udata = Udata;

   % Local copies, udata{x}(4) is just too complicated...
   data            = udata{1};
   tmp        = udata{2};
   lolim = tmp(1);
   hilim = tmp(2);
   l1              = udata{3};
   l2              = udata{4};
   l3              = udata{5};
   tmp = udata{6};
   xi1 = tmp(1);
   xi2 = tmp(2);
   xi3 = tmp(3);
   % udata{7} -> axes of subplots, get from subplot();
   % udata{8} -> function handle to this function, not used
end

if (isempty(udata{7}))
   % old_axis must have be valid axis.  Create one since
   % no subplots exist yet in this figure

   subplot(2,2,1)
end

old_axis = gca;

himg = subplot(2,2,1);
imagesc(l1, l2, data(:,:,xi3), [lolim hilim]); axis image;
set(gca,'YDir','Normal');
xlabel('X'); ylabel('Y'); title(sprintf('HSI5: Z = %d', xi3));
crosshair(l1(1), l1(xi1), l1(end), l2(1), l2(xi2), l2(end));
% colorbar;

hcol = subplot(2,2,2);
imagesc(l3, l2, squeeze(data(:,xi1,:)), [lolim hilim]); axis image;
set(gca, 'YDir', 'Normal', 'YAxisLocation', 'Right');
xlabel('Z'); ylabel('Y'); title(sprintf('X = %d', xi1));
crosshair(l3(1), l3(xi3), l3(end), l2(1), l2(xi2), l2(end));

hrow = subplot(2,2,3);
imagesc(l1, l3, squeeze(data(xi2,:,:))', [lolim hilim]); axis image;
set(gca,'YDir','Normal');
xlabel('X'); ylabel('Z'); title(sprintf('Y = %d',xi2));
crosshair(l1(1), l1(xi1), l1(end), l3(1), l3(xi3), l3(end));

hpix = subplot(2,2,4);
plot(l3, squeeze(data(xi2,xi1,:)));
axis([l3(1) l3(end) lolim hilim]);
xlabel('Z'); title(sprintf('[X,Y] = [%d,%d]', xi1, xi2));

% Copy axis back for caller
udata{7} = [ himg, hcol, hrow, hpix ];

% Restore original 'current' axis
axes(old_axis);
%drawnow;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the plot from inside the callback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[udata] = updateScreen(Udata)

% Initialize with previous values

udata = Udata;

% Get current point

pt = get(gca, 'CurrentPoint');

hin = pt(1,1);
vin = pt(1,2);
% disp([ gca hin vin ])

% Force point to remain within axis limits

hlim = get(gca, 'XLim');
hin = max(hlim(1), min(hlim(2), hin));

vlim = get(gca, 'YLim');
vin = max(vlim(1), min(vlim(2), vin));

% Convert from axis limits to array index

h1   = gca;
hsub = udata{7};

if (isempty(hsub))
   error('Empty subplot list');
   return;
end

l1 = udata{3};   n1 = length(l1);
l2 = udata{4};   n2 = length(l2);
l3 = udata{5};   n3 = length(l3);

tmp = udata{6};
x1 = tmp(1);
x2 = tmp(2);
x3 = tmp(3);

switch (h1)
   case hsub(1)
      x1p = findindex(l1, hin);
      x2p = findindex(l2, vin);
      x3p = x3;
   case hsub(2)
      x1p = x1;
      x2p = findindex(l2, vin);
      x3p = findindex(l3, hin);
   case hsub(3)
      x1p = findindex(l1, hin);
      x2p = x2;
      x3p = findindex(l3, vin);
   case hsub(4)
      % Motion in Z only
      x1p = x1;
      x2p = x2;
      x3p = findindex(l3, hin);
   otherwise
      error('Not an image axis');
end

% Force index back into legal range, just in case

x1 = max(1, min(n1, x1p));
x2 = max(1, min(n2, x2p));
x3 = max(1, min(n3, x3p));

% Update display using new coordinates

udata{6} = [ x1 x2 x3 ];

udata = feval(udata{8}, udata);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw cross-hairs across plot marking the other image planes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function crosshair(u1,ui,u2,v1,vi,v2)

hold on;

plot([u1,u2],[vi,vi], 'w');
plot([ui,ui],[v1,v2], 'w');

hold off;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn axis coordinates into matrix indices
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
