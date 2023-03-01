% getPhiTotalN - return PhiTotalN for a specified wavelength in a
%              vector (amplitude and phase) that is a 1-to-1 match
%              with the MeasList
%
%   v = getPhiTotalN( pmiModel,idxLambda )
%
%      pmiModel     The model describing the measurement list and
%                   modulation frequencies.
%
%      idxLambda    Specifies the wavelength
%

%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/04/27 16:23:41 $
%
%  $Revision: 1.1 $
%
%  $Log: getPhiTotalN.m,v $
%  Revision 1.1  2001/04/27 16:23:41  dboas
%  Initial version
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = getPhiTotalN( pmiModel, idxLambda )

  list0 = find( pmiModel.MeasList(:,4)==idxLambda );
  
  nMeas = length( list0 );
  v = zeros(nMeas,1);
  
  
  phi = pmiModel.P(idxLambda).PhiTotalN;
  
  nFreq = length( pmiModel.ModFreq );
  
  
  idxMeas = 1;
  for idxFreq = 1:nFreq
    
    list = find(pmiModel.MeasList(list0,3)==idxFreq );
    nMeasF = length(list);
    
    if pmiModel.ModFreq(idxFreq)==0
      v(list) = phi(idxMeas:idxMeas+nMeasF-1);
      idxMeas = idxMeas + nMeasF;
    else
      v(list) = phi(idxMeas:idxMeas+nMeasF-1) + ...
		sqrt(-1) * phi(idxMeas+nMeasF:idxMeas+2*nMeasF-1);
      idxMeas = idxMeas + 2*nMeasF;
    end
    
  end
  
    