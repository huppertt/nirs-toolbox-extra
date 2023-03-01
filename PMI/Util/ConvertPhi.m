% ConvertPhi - Convert the P.Phi vector into a complex
%              vector (amplitude and phase) that is a 1-to-1 match
%              with the MeasList
%
%   v = ConvertPhi( pmiModel, phi, idxLambda )
%
%      pmiModel     The model describing the measurement list and
%                   modulation frequencies.
%
%      phi          The fluence, either P(idxLambda).PhiTotalN,
%                   P(idxLambda).PhiPhiScat, etc.
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/04/27 15:13:38 $
%
%  $Revision: 1.0 $
%
%  $Log: ConvertPhi,v $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = ConvertPhi( pmiModel, phi, idxLambda )

  nMeas = length( find(pmiModel.MeasList(:,4)==idxLambda) );
  nFreq = length( pmiModel.ModFreq );
  
  v = zeros(nMeas,1);

  idxMeas = 1;
  for idxFreq = 1:nFreq
    
    list = find(pmiModel.MeasList(:,3)==idxFreq & ...
		pmiModel.MeasList(:,4)==idxLambda );
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
  