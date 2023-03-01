%multImage      Interactive display of a number of PMI images.
%
%   idxCurrFig = mutlImage(ds, idxCurrFig, xEst)
%
%   ds          The DPDW Imaging data structure containing the
%               reconstruction to be imaged.
%
%   idxCurrFig  OPTIONAL: The number of the figure in which to display the
%               images.  It is also return incremented by the number of
%               figures displayed.
%
%   xEst        OPTIONAL: The reconstructions to display if they are not in
%               ds.Recon.xExt.
%
%   multImage displays the reconstructions present in the PMI imaging
%   data structure or the separately supplied 
%
%   Calls: getvisinfo, mcontour, islice.
%
%   Bugs: only slice type images are implemented, since multple wavelengths
%   are indexed by the same dimension as multiple reconstructions.

% Copyright (C) 2002, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: multImage.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 19:29:38  rjg
%  Initial Revision for PMI 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function idxCurrFig = multImage(ds, idxCurrFig, xEst)
if nargin < 2
    idxCurrFig = gcf;
end
if nargin > 2
    ds.Recon.xEst = xEst;
end



figure(idxCurrFig)
clf
orient landscape

%%
%%  Reconstruct image of the volume
%%
switch ds.Visualize.CMap
case 'grey'
    colormap(flipud(gray(256)));
case 'bgyor'
    colormap(jet(256));
end

switch ds.Visualize.Type
case 'contour'
    disp('Not yet implemented');
case 'image'
    switch ds.Visualize.CRange
    case 'fixed'
        islice(ds.Recon.xEst, ds.Inv.CompVol, ...
            ds.Visualize.PlaneIndices, ...
            ds.Visualize.VisPlane, ...
            ds.Visualize.LayoutVector, ds.Visualize.CRangeVal);
                    
    case 'auto'
        islice(ds.Recon.xEst, ds.Inv.CompVol, ...
            ds.Visualize.PlaneIndices, ...
            ds.Visualize.VisPlane, ds.Visualize.LayoutVector);
    end
end

idxCurrFig = idxCurrFig + 1;
