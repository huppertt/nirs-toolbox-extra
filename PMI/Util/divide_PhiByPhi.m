% divide_PhiByPhi - return the complex ratio of phi1 by phi2 after
%                   appropriately recombining phi1 and phi2 into
%                   complex numbers if necessary.
%
%   v = divide_PhiByPhi( pmiModel, phi1, phi2, idxLambda );
%
%      pmiModel     The model describing the measurement list and
%                   modulation frequencies.
%
%      phi1, phi2   The data corresponding to the measurement list.
%
%     idxLambda     Indicates which lambda to use.
%
%   Calls: none.
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/08/16 15:13:38 $
%
%  $Revision: 1.0 $
%
%  $Log: divide_PhiByPhi,v $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = divide_PhiByPhi( pmiModel, phi1, phi2, idxLambda )
  
  nMeas = length( find(pmiModel.MeasList(:,4)==idxLambda) );
  v = zeros(nMeas,1);
  
  for idxMeas = 1:nMeas

    if pmiModel.ModFreq == 0
      v(idxMeas) = phi1(idxMeas) / phi2(idxMeas);
    else
      v(idxMeas) = (phi1(idxMeas) + i * phi1(idxMeas+nMeas)) / ...
	  (phi2(idxMeas) + i * phi2(idxMeas+nMeas));
    end
    
  end
  