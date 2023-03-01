%MSLICE         Show multiple slices of a 3D volume.
%
%   mslice(X, CompVol, idxSlices, VisPlane, dimPlot, CAxis)
%
%   X           The 3D volume to show.
%
%   CompVol     A structure defining the computational volume.  This
%               structure should have the members: Type, X, Y and Z.  Type
%               should be uniform specifying a uniform sampling volume of
%               voxels. X, Y and Z are vectors specifying the centers of the
%               voxels.
%
%   idxSlices   OPTIONAL: The slices to show.
%
%   VisPlane    OPTIONAL: The dimension in which to slice up the volume.
%               The visualization plane will be constant in this
%               dimension.
%
%   dimPlot     OPTIONAL: The subplot dimensions (default: 2x3).
%
%   CAxis       OPTIONAL: The color axis range for all of the plots.  The
%               default is the range of each slice independently.
%
%   MSLICE
%
%   Calls: none.
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:48 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: mslice.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:48  dboas
%  initial
%
%  Revision 2.4  1999/11/15 18:22:29  rjg
%  Moved to PMI.
%
%  Revision 2.3  1999/02/05 21:06:47  rjg
%  Correc axes permutation problem and handling of color range
%
%  Revision 2.2  1998/08/20 19:39:27  rjg
%  Fixed help to reflect VisPlane parameter.
%
%  Revision 2.1  1998/08/07 21:28:49  rjg
%  Added the ability to slice in all three directions.
%
%  Revision 2.0  1998/08/05 18:30:22  rjg
%  Handles CompVol structure as part of the SlabImage data structure.
%
%  Revision 1.1  1998/04/29 16:00:03  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mslice(X, CompVol, idxSlices, VisPlane, dimPlot, CAxis)

nX = length(CompVol.X);
nY = length(CompVol.Y);
nZ = length(CompVol.Z);

%%
%%  Reshape the data in vector format, this won't affect an array
%%  already in 3D.
%%
X = reshape(X, nY, nX, nZ);

%%
%%  If the slice indices is not already given subsample so that there are 6
%%  images uniform located in depth.
%%
if nargin < 5
    nRows = 2;
    nCols = 3;
    if nargin < 4
        VisPlane = 'Z';
        if nargin < 3
            idxSlices = [1:nZ/6:nZ];
            idxSlices = round(idxSlices);
        end
    end
else
    nRows = dimPlot(1);
    nCols = dimPlot(2);
end

switch VisPlane
case 'X'
    ColValues = CompVol.Y;
    RowValues = CompVol.Z;
    SliceMax = size(X, 1);
case 'Y'
    ColValues = CompVol.X;
    RowValues = CompVol.Z;
    SliceMax = size(X, 2);
case 'Z'
    ColValues = CompVol.X;
    RowValues = CompVol.Y;
    SliceMax = size(X, 3);
end


cmin = min(X(:));
cmax = max(X(:));
if nargin > 5
    if CAxis(1) == CAxis(2)
        CAxis = [cmin cmax];
    end
end

for idxPlot = 1:length(idxSlices)
    subplot(nRows, nCols, idxPlot);

    switch VisPlane
    case 'X'
        imagesc(ColValues, RowValues, squeeze(X(:,idxSlices(idxPlot),:))');
        title(['X = ' num2str(CompVol.X(idxSlices(idxPlot)))]);
    case 'Y'
        imagesc(ColValues, RowValues, squeeze(X(idxSlices(idxPlot),:,:))');
        title(['Y = ' num2str(CompVol.Y(idxSlices(idxPlot)))]);
    case 'Z'
        imagesc(ColValues, RowValues, X(:,:, idxSlices(idxPlot)));
        title(['Z = ' num2str(CompVol.Z(idxSlices(idxPlot)))]);
    end
    if nargin > 5,
        caxis(CAxis);
    end
    grid on
    set(gca, 'ydir', 'normal');
    colorbar('vert');

    if idxPlot == SliceMax
        break
    end
end

