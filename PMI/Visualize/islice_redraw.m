%ISLICE_REDRAW  Redraw function for ISLICE
%
%   islice_redraw
%
%
%   ISLICE_REDRAW
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
%  $Log: islice_redraw.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:48  dboas
%  initial
%
%  Revision 1.0  1999/11/15 18:17:25  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function islice_redraw

IDS = get(gcf, 'UserData');

%%
%%  Reshape the data in vector format, this won't affect an array
%%  already in 3D.
%%
iCurr = get(IDS.Slider, 'value') * IDS.nDataCols;
if iCurr < 1
    iCurr = 1;
end
if iCurr > IDS.nDataCols
    iCurr = IDS.nDataCols;
end
iCurr = floor(iCurr);
set(IDS.Text, 'string', int2str(iCurr));
X = reshape(IDS.objfcts(:, iCurr), IDS.nY, IDS.nX, IDS.nZ);

if IDS.flgCConstant
    cmin = min(X(:));
    cmax = max(X(:));
    IDS.CAxis = [cmin cmax];
end

for idxPlot = 1:length(IDS.idxSlices)
    subplot(IDS.nRows, IDS.nCols, idxPlot);

    switch IDS.VisPlane
    case 'X'
        imagesc(IDS.ColValues, IDS.RowValues, ...
            squeeze(X(:,IDS.idxSlices(idxPlot),:))');
        set(gca, 'ydir', 'normal');
        title(['X = ' num2str(IDS.CompVol.X(IDS.idxSlices(idxPlot)))]);
    
    case 'Y'
        imagesc(IDS.ColValues, IDS.RowValues, ...
            squeeze(X(IDS.idxSlices(idxPlot),:,:))');
        set(gca, 'ydir', 'normal');
        title(['Y = ' num2str(IDS.CompVol.Y(IDS.idxSlices(idxPlot)))]);
    
    case 'Z'
        imagesc(IDS.ColValues, IDS.RowValues, X(:,:, ...
            IDS.idxSlices(idxPlot)));
        set(gca, 'ydir', 'normal');
        title(['Z = ' num2str(IDS.CompVol.Z(IDS.idxSlices(idxPlot)))]);
    end
    if ~IDS.flgCAuto
        caxis(IDS.CAxis);
    end
    grid on
    colorbar('vert');

    if idxPlot == IDS.SliceMax
        break
    end
end

