%MCONTOUR       Show multiple contours of a 3D volume.
%
%   mcontour(X, CompVol, idxSlices, VisPlane, dimPlot, nCLines)
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
%   nCLines     OPTIONAL: The values of the contour lines to display.
%
%   MCONTOUR
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
%  $Log: mcontour.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:48  dboas
%  initial
%
%  Revision 2.2  1999/11/15 18:18:24  rjg
%  Moved to PMI.
%
%  Revision 2.1  1999/02/05 21:04:19  rjg
%  Added visualization plane selection.
%
%  Revision 2.0  1998/08/05 18:27:03  rjg
%  Handles CompVol structure as part of the SlabImage data structure.
%
%  Revision 1.1  1998/04/29 16:02:51  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mcontour(X, CompVol, idxSlices, VisPlane, dimPlot, nCLines)

nX = length(CompVol.X);
nY = length(CompVol.Y);
nZ = length(CompVol.Z);

%%
%%  Reshape the data in vector format, this won't affect a array
%%  already in 3D.
%%
X = reshape(X, nX, nY, nZ);

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

for idxPlot = 1:length(idxSlices)
    subplot(nRows, nCols, idxPlot);

    switch VisPlane
    
    case 'X'
        if nargin > 5
            cs = contour(ColValues, RowValues, ...
                squeeze(X(:,idxSlices(idxPlot),:))', nCLines);
        else
            cs = contour(ColValues, RowValues, ...
                squeeze(X(:,idxSlices(idxPlot),:))');
        end
        title(['X = ' num2str(CompVol.X(idxSlices(idxPlot)))]);
    
    case 'Y'
        if nargin > 5
                    cs = contour(ColValues, RowValues, ...
                squeeze(X(idxSlices(idxPlot),:,:))', nCLines);
        else
            cs = contour(ColValues, RowValues, ...
                squeeze(X(idxSlices(idxPlot),:,:))');
        end
        title(['Y = ' num2str(CompVol.Y(idxSlices(idxPlot)))]);
    
    case 'Z'
        if nargin > 5
            cs = contour(ColValues, RowValues, X(:,:, idxSlices(idxPlot)), ...
                nCLines);
        else
            cs = contour(ColValues, RowValues, X(:,:, idxSlices(idxPlot)));
        end
        title(['Z = ' num2str(CompVol.Z(idxSlices(idxPlot)))]);
    end

    clabel(cs)
    grid on
end

