% PLOTDATA  2D plots of the data based on source and detector location,
%           one subplot per detector.
%
% h = plotData(SD, MeasList, Phi, srcFlag, cmap);
%
% Inputs: 
%       SD - SD structure, used for SrcPos, DetPos
% MeasList - data elements to be plotted
%      Phi - The data to plot.  Should match MeasList.
%  srcFlag - OPTIONAL: If non-zero, plot per-source instead of per-detector.
%     cmap - OPTIONAL: colormap to use
%
% Outputs:
%    h   - Figure handle
%
% BUGS: Multiple wavelengths/frequencies/times are not properly 
%       supported correctly.  If a src-det pair occurs in the measurement 
%       list more than once, the initial data will be over-written.

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

function[h] = plotData(SD, MeasList, PhiMeas, doSrc, cmap);

optodeColor = 'magenta';

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('doSrc','var') | isempty(doSrc))
  doSrc = 0;
end

h  = figure;

if (~exist('cmap','var') | isempty(cmap))
   if (exist('cm256'))
     % Use my nice colormap, if available
     colormap(cm256);
   end
else
   if (size(cmap,2) ~= 3)
      error('cmap must be an Nx3 matrix');
   else
      colormap(cmap)     
   end
end

PhiMeas = PhiMeas(:);

% If data sizes don't match, assume that only one wavelength (the first
%  one) is present.

if (size(MeasList,1) > length(PhiMeas))
   MeasList = MeasList( find(MeasList(:,4)==1), :);
end

% Compute image boundaries

x1 = min( min(SD.SrcPos(MeasList(:,1),1)), ...
	  min(SD.DetPos(MeasList(:,2),1)) );
x2 = max( max(SD.SrcPos(MeasList(:,1),1)), ...
	  max(SD.DetPos(MeasList(:,2),1)) );

y1 = min( min(SD.SrcPos(MeasList(:,1),2)), ...
	  min(SD.DetPos(MeasList(:,2),2)) );
y2 = max( max(SD.SrcPos(MeasList(:,1),2)), ...
	  max(SD.DetPos(MeasList(:,2),2)) );

z1 = min(PhiMeas);
z2 = max(PhiMeas);

if (z1 >= z2)
   warning('z1 >= z2');

   % Data vector is constant, broaden range a bit
   z1 = z1 - 1;
   z2 = z2 + 1;
end

% Pad out boundaries

x1 = floor(x1 - 1); x2 =  ceil(x2 + 1);
y1 = floor(y1 - 1); y2 =  ceil(y2 + 1);

cm = colormap;
nc = size(cm,1);

if (doSrc)
   dlst = unique(MeasList(:,2));

   % Number of subplots needed
   % nf = ceil(sqrt(size(SD.DetPos,1)));
   [M,N] = plotDimensions(length(dlst), 1, 0);

   for iD = 1:length(dlst)
      idxD = dlst(iD);
      
      subplot(M,N,iD);
      axis([x1 x2 y1 y2]);
	 
      ml = find(MeasList(:,2) == idxD);

      if (~isempty(ml))
	 % Plot each source measurement of this detector as a patch
	 
	 for idxM = 1:length(ml)
	    si = MeasList(ml(idxM),1);
	    sx = SD.SrcPos(si,1);
	    sy = SD.SrcPos(si,2);
	    
	    ci = round((nc-1) * (PhiMeas(ml(idxM)) - z1) / (z2 - z1));
	    
	    sz = cm(ci+1,:);
	    
	    patch( [ sx-0.5, sx+0.5, sx+0.5, sx-0.5, sx-0.5 ], ...
		   [ sy-0.5, sy-0.5, sy+0.5, sy+0.5, sy-0.5 ], sz);
	 end
	 
	 dx = SD.DetPos(idxD,1);
	 dy = SD.DetPos(idxD,2);
	 
	 line([ dx-0.5 dx+0.5 ], [ dy-0.5 dy+0.5 ], 'Color', optodeColor);
	 line([ dx-0.5 dx+0.5 ], [ dy+0.5 dy-0.5 ], 'Color', optodeColor);
	 line([ dx-0.5 dx+0.5 ], [ dy-0.5 dy+0.5 ], 'Color', optodeColor);
	 line([ dx-0.5 dx+0.5 ], [ dy+0.5 dy-0.5 ], 'Color', optodeColor);
      end
      
      ylabel([ 'Detector ' num2str(idxD) ]);
      
      if (isempty(ml))
	 axis off
      end
   end
else
   slst = unique(MeasList(:,1));
   
   % Number of subplots needed
   %nf = ceil(sqrt(size(SD.SrcPos,1)));
   [M,N] = plotDimensions(length(slst), 1, 0);
   
   for iS = 1:length(slst)
      idxS = slst(iS);
      
      subplot(M,N,iS);
      axis([x1 x2 y1 y2]);
   
      ml = find(MeasList(:,1) == idxS);
   
      if (~isempty(ml))
	 % Plot each source measurement of this detector as a patch
	 
	 for idxM = 1:length(ml)
	    di = MeasList(ml(idxM),2);
	    dx = SD.DetPos(di,1);
	    dy = SD.DetPos(di,2);
	    
	    ci = round((nc-1) * (PhiMeas(ml(idxM)) - z1) / (z2 - z1));
	    
	    dz = cm(ci+1,:);
	    
	    patch( [ dx-0.5, dx+0.5, dx+0.5, dx-0.5, dx-0.5 ], ...
		   [ dy-0.5, dy-0.5, dy+0.5, dy+0.5, dy-0.5 ], dz);
	 end
      
	 sx = SD.SrcPos(idxS,1);
	 sy = SD.SrcPos(idxS,2);
   
	 line([ sx-0.5 sx+0.5 ], [ sy-0.5 sy+0.5 ], 'Color', optodeColor);
	 line([ sx-0.5 sx+0.5 ], [ sy+0.5 sy-0.5 ], 'Color', optodeColor);
      end   
      
      ylabel([ 'Source ' num2str(idxS) ]);
      
      if (isempty(ml))
	 axis off
      end
   end
end

return;
