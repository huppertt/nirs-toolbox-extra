%ISLICE         Interactive display of multiple slices of a 3D volume.
%
%   islice(objfcts, CompVol, idxSlices, VisPlane, dimPlot, CAxis)
%
%   objfcts     The 3D volumes to show.
%
%   CompVol     A structure defining the computational volume.  This
%               structure should have the members: Type, X, Y and Z.  Type
%               should be uniform specifying a uniform sampling volume of
%               voxels. X, Y and Z are vectors specifying the centers of the
%               voxels.
%
%   idxSlices   OPTIONAL: The slices to show (default [11:-1:1]).
%
%   VisPlane    OPTIONAL: The dimension in which to slice up the volume.
%               The visualization plane will be constant in this
%               dimension.
%
%   dimPlot     OPTIONAL: The subplot dimensions (default: 2x3).
%
%   CAxis       OPTIONAL: The color axis range for all of the plots.  The
%               default is the total range of each object function.  A
%               single 0 will auto range each slice independently.
%
%   ISLICE
%
%   Calls: islice_redraw.
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
%  $Log: islice.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:48  dboas
%  initial
%
%  Revision 1.1  1999/11/15 18:16:40  rjg
%  Fixed RCS author field
%
%  Revision 1.0  1999/11/15 18:15:38  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function islice(objfcts, CompVol, idxSlices, VisPlane, dimPlot, CAxis)

IDS.objfcts = objfcts;
IDS.CompVol = CompVol;
IDS.nX = length(CompVol.X);
IDS.nY = length(CompVol.Y);
IDS.nZ = length(CompVol.Z);

%%
%%  If the slice indices is not already given subsample so that there are 6
%%  images uniform located in depth.
%%
if nargin < 5
    IDS.nRows = 3;
    IDS.nCols = 3;
else
    IDS.nRows = dimPlot(1);
    IDS.nCols = dimPlot(2);
end
if nargin < 4
    IDS.VisPlane = 'Z';
else
    IDS.VisPlane = VisPlane;
end
if nargin < 3
    IDS.idxSlices = [9:-1:1];
else
    IDS.idxSlices = idxSlices;
end

switch IDS.VisPlane
case 'X'
    IDS.ColValues = IDS.CompVol.Y;
    IDS.RowValues = IDS.CompVol.Z;
    IDS.SliceMax = length(IDS.CompVol.X);
case 'Y'
    IDS.ColValues = IDS.CompVol.X;
    IDS.RowValues = IDS.CompVol.Z;
    IDS.SliceMax = length(IDS.CompVol.Y);
case 'Z'
    IDS.ColValues = IDS.CompVol.X;
    IDS.RowValues = IDS.CompVol.Y;
    IDS.SliceMax = length(IDS.CompVol.Z);
end

IDS.flgCConstant = 0;
IDS.flgCAuto = 0;
if nargin > 5
    if length(CAxis) == 1 & CAxis == 0
        IDS.flgCAuto = 1;
    else
        IDS.CAxis = CAxis;
    end
else
    IDS.flgCConstant = 1;
end


%%
%%  Setup the slider control
%%
IDS.nDataCols = size(objfcts,2);
if IDS.nDataCols > 10
    sstep = [1/IDS.nDataCols 10/IDS.nDataCols];
else
    sstep = [1/IDS.nDataCols 1];
end

IDS.Slider = uicontrol(	'Units', 'Normalized', ...
    'Position', [0.05 0.1 0.025 0.8], ...
	'Style','slider', ...
    'SliderStep', sstep, ...
    'callback', 'islice_redraw');
IDS.Text= uicontrol('Units','Normalized', ...
	'Position',[0.05 0.9 0.025 0.028], ...
	'Style','text');

set(gcf, 'UserData', IDS);
islice_redraw
