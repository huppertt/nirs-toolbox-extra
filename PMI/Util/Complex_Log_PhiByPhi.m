% divide_PhiByPhi - return the log of the complex ratio of phi1 by phi2 after
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
%  $Date: 2000/11/08 18:08:44 $
%
%  $Revision: 1.2 $
%
%  $Log: Complex_Log_PhiByPhi.m,v $
%  Revision 1.2  2000/11/08 18:08:44  dboas
%  Fixed this routine to properly handle data sets with multiple
%  modulation frequencies.
%
%  Revision 1.1  2000/09/28 20:06:40  dboas
%  Initial version.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = Complex_Log_PhiByPhi( pmiModel, phi1, phi2, idxLambda )
 
  nFreq = length(pmiModel.ModFreq);
  nLambda = length(pmiModel.Lambda);

  nData = size(pmiModel.MeasList,1)/(nLambda*nFreq) * ...
    (nFreq + sum(pmiModel.ModFreq > 0));

  v = zeros(nData,1);

  idxRow = 0;
  for idxFreq = 1:nFreq

    FreqList = find(pmiModel.MeasList(:,3)==idxFreq);
 
    nMeas = length( find(pmiModel.MeasList(FreqList,4)==idxLambda) );

    for idxMeas = 1:nMeas

      if pmiModel.ModFreq(idxFreq) == 0
        v(idxRow+idxMeas) = log(phi1(idxRow+idxMeas) / phi2(idxRow+idxMeas));
      else
        foo = log( (phi1(idxRow+idxMeas) + i * phi1(idxRow+idxMeas+nMeas)) / ...
	    (phi2(idxRow+idxMeas) + i * phi2(idxRow+idxMeas+nMeas)) );
        v(idxRow+idxMeas) = real(foo);
        v(idxRow+idxMeas+nMeas) = imag(foo);
      end

    end

    if pmiModel.ModFreq(idxFreq) == 0
      idxRow = idxRow + nMeas;
    else
      idxRow = idxRow + 2*nMeas;
    end

  end    