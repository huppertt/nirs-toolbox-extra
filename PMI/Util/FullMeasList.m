%FullMeasList   Create a full measurement list in standard form.
%
%   MeasList = FullMeasList(nSrc, nDet, nFreq, nLambda)
%
%   MeasList    The MeasList matrix.
%
%   nSrc, nDet, nFreq, nLambda  The number of each measurement parameter for the
%                               measurement list.
%
%   FullMeasList returns a PMI Measurement List field specfying that all
%   detectors are cycled through first, then all sources then all frequencies
%   then finally all wavelengths.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: FullMeasList.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.1  1999/11/10 18:50:32  rjg
%  Swapped the positions of the source and detector indices in the MeasList
%  result and in the calling format.
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MeasList = FullMeasList(nSrc, nDet, nFreq, nLambda)
%%
%%  Create the full MeasList matrix
%%
DOrder = repmat([1:nDet]', nSrc*nFreq*nLambda, 1);
tmp = repmat([1:nSrc], nDet, 1);
tmp = tmp(:);
SOrder = repmat(tmp, nFreq*nLambda, 1);
tmp = repmat([1:nFreq], nDet*nSrc, 1);
tmp = tmp(:);
FOrder = repmat(tmp, nLambda, 1);
tmp = repmat([1:nLambda], nDet*nSrc*nFreq, 1);
LOrder = tmp(:);

MeasList = [SOrder DOrder FOrder LOrder];
