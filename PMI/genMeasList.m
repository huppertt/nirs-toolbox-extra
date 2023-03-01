%GenMeasList    Generate a subset measurement list from the rawdata measurement 
%               list and return it.
%
%   MeasList = genMeasList(rawData, method, var1...);
%
%   rawData     This is the raw data time series structure.
%
%   method      This defines the method used to create a subset
%               measurement list.  Presently it can equal
%               'Nearest' - Choose the nearest neighbor source-detector
%               	pairs such that the separation is less than var1.
%					 'SNR' - Choose the threshold SNR.
%						var1 - threshold SNR
%						var2 - data frame containing SNR
%
%   var1        A variable whose use depends on METHOD.
%
%
%   Calls: none.
%
%   Bugs: none.

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
%  $Date: 2001/04/19 20:27:56 $
%
%  $Revision: 1.5 $
%
%  $Log: genMeasList.m,v $
%  Revision 1.5  2001/04/19 20:27:56  dboas
%  SNR was fixed to find the correspondence between a given MeasList
%  and pmi.data.MeasList
%
%  Revision 1.4  2000/07/27 14:56:04  dboas
%  Removed a bunch of old commented lines.
%
%  Revision 1.3  2000/06/27 19:30:54  dboas
%  Whoops, now SNR is fixed.
%
%  Revision 1.2  2000/06/27 19:25:35  dboas
%  Fixed the 'SNR' selection to make it work.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.3  2000/01/10 00:14:14  dboas
%  Storing the source and detector lists for use by other functions
%
%  Revision 3.2  1999/12/03 13:53:10  dboas
%  Fixed a small typo error from previous submission
%
%  Revision 3.1  1999/11/16 23:38:25  dab
%  Initial revision.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MeasList] = genMeasList( method, var0, var1, var2, var3 )

switch lower(method)
 
 case 'nearest'
   
  nMeasList = 0;
  for ii=1:length(var0.MeasList)
    rs = squeeze(var0.Src.Pos(var0.MeasList(ii,1),:));
    rd = squeeze(var0.Det.Pos(var0.MeasList(ii,2),:));
    r = ( (rs-rd) * (rs-rd)' )^0.5;
    
    if r<var1
      nMeasList = nMeasList + 1;
      index(nMeasList) = ii;
    end
  end
  
  if nMeasList>0
    MeasList = var0.MeasList(index,:);
  end
  
 case 'snr'
  nMeasList = 0;
  for ii=1:length(var0)
    jj = find(var0(ii,1)==var1.MeasList(:,1) & ...
	      var0(ii,2)==var1.MeasList(:,2) & ...
	      var0(ii,3)==var1.MeasList(:,3) & ...
	      var0(ii,4)==var1.MeasList(:,4) );
    SNR = abs(var1.Raw(jj,var3) / var1.Raw_std(jj,var3));
    if SNR>var2
      nMeasList = nMeasList +1;
      index(nMeasList) = ii;
    end
  end

  if nMeasList>0
    MeasList = var0(index,:);
  end
  
 case 'all'
  nFreq = size(var0.ModFreq,2);
  nLambda = size(var0.Lambda,2);
  [pSrc nSrc] = getOptodePos(var0.Src);
  [pDet nDet] = getOptodePos(var0.Det);    
  MeasList = FullMeasList(nSrc, nDet, nFreq, nLambda);

 otherwise
  disp('Unknown method for generating a measurement list.');
  
end

   