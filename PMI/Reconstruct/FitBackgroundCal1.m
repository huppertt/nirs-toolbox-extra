%function [delta_mu, lsq, CCs, CCd, PhiInc]=FitBackgroundCal1(Mu,pmiModel,idxLambda, lambda)
%
% This is a sub-function to FitBackgroundCal and is used to get
% the Newton-Raphson update
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/04/19 19:48:40 $
%
%  $Revision: 1.5 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta_mu, lsq, CCs, CCd, PhiInc]=FitBackgroundCal1(Mu,pmiModel,idxLambda, lambda)

  delta_musp = 0.01;
  delta_mua = 0.0001;
  
  idxThisLambda = find(idxLambda == pmiModel.MeasList(:,4));

  if pmiModel.ModFreq == 0
    lnPhiM = log(pmiModel.P(idxLambda).PhiTotalN);
  else
    nMeas = size(find(pmiModel.MeasList(idxThisLambda,4)==idxLambda),1);
    clear foo;
    foo = pmiModel.P(idxLambda).PhiTotalN;
    lnPhiM(1:nMeas,1) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnPhiM(nMeas+1:nMeas*2,1) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
  end

  nSrcs = 0;
%  srcIdx = zeros(length(pmiModel.Src.Pos),1);
  for idx = 1:length(pmiModel.Src.Pos)
    clear foo;
    foo = find( idx==pmiModel.MeasList(idxThisLambda,1) );
    if ~isempty( foo )
      nSrcs = nSrcs + 1;
      srcIdx(idx) = nSrcs;
    end
  end
  nDets = 0;
  for idx = 1:length(pmiModel.Det.Pos)
    clear foo;
    foo = find( idx==pmiModel.MeasList(idxThisLambda,2) );
    if ~isempty( foo )
      nDets = nDets + 1;
      detIdx(idx) = nDets;
    end
  end
  nCC  = nSrcs + nDets - 1;

  MeasList = pmiModel.MeasList(idxThisLambda,:);
  nMeas = size(MeasList,1);

  %determine the calibration jacobian associated with the source
  %and detector coupling coefficients
  if pmiModel.ModFreq == 0
    C = zeros( nMeas, nCC);
    for i=1:nMeas
      if srcIdx(MeasList(i,1)) ~= 1
        C(i,srcIdx(MeasList(i,1))-1) = 1;  %This is the source (
					  %minus 1 because we know
					  %the first source
                                          %coupling coefficient)
      end
      C(i,detIdx(MeasList(i,2))+nSrcs-1) = 1;
    end
  else
    C = zeros( 2*nMeas, 2*nCC);
    for i=1:nMeas
      if srcIdx(MeasList(i,1)) ~= 1
        C(i,srcIdx(MeasList(i,1))-1) = 1;  %This is the source (
					  %minus 1 because we know
					  %the first source
                                          %coupling coefficient)
        C(i+nMeas,srcIdx(MeasList(i,1))+nCC-1) = 1;
      end
      C(i,detIdx(MeasList(i,2))+nSrcs-1) = 1;
      C(i+nMeas,detIdx(MeasList(i,2))+nCC+nSrcs-1) = 1;
    end
  end
  
  %Evaluate at Mu_o
  pmiModel.Mu_sp(idxLambda) = Mu(1);
  pmiModel.Mu_a(idxLambda) = Mu(2);
  FwdModel=genMeasData(pmiModel);
  PhiInc = FwdModel.P(idxLambda).PhiInc;
  if pmiModel.ModFreq == 0
    lnPhiTo = log(FwdModel.P(idxLambda).PhiTotal);
    lnCC = -(C'*C) \ (C' * (lnPhiTo - lnPhiM));
    lnPhiT = lnPhiTo + C*lnCC;
    lnCCo = lnCC;
  else
    clear foo;
    foo = FwdModel.P(idxLambda).PhiTotal;
    lnPhiTo(1:nMeas,1) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnPhiTo(nMeas+1:nMeas*2,1) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
    lnCC = -(C'*C) \ (C' * (lnPhiTo - lnPhiM));
    lnPhiT = lnPhiTo + C*lnCC;
    lnCCo = lnCC;
  end

  %Evaluate the derivative versus mu_sp
  pmiModel.Mu_sp(idxLambda) = Mu(1)+delta_musp;
  pmiModel.Mu_a(idxLambda) = Mu(2);
  FwdModel=genMeasData(pmiModel);
  if pmiModel.ModFreq == 0
    lnPhiTo = log(FwdModel.P(idxLambda).PhiTotal);
    lnCC = -(C'*C) \ (C' * (lnPhiTo - lnPhiM));
    dlnPhiT_dmusp = (lnPhiTo+C*lnCC-lnPhiT) / delta_musp;
  else
    clear foo;
    foo = FwdModel.P(idxLambda).PhiTotal;
    lnPhiTo(1:nMeas,1) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnPhiTo(nMeas+1:nMeas*2,1) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
    lnCC = -(C'*C) \ (C' * (lnPhiTo - lnPhiM));
    dlnPhiT_dmusp = (lnPhiTo+C*lnCC-lnPhiT) / delta_musp;
  end
  
  %Evaluate the derivative versus mu_a
  pmiModel.Mu_sp(idxLambda) = Mu(1);
  pmiModel.Mu_a(idxLambda) = Mu(2) + delta_mua;
  FwdModel=genMeasData(pmiModel);
  if pmiModel.ModFreq == 0
    lnPhiTo = log(FwdModel.P(idxLambda).PhiTotal);
    lnCC = -(C'*C) \ (C' * (lnPhiTo - lnPhiM));
    dlnPhiT_dmua = (lnPhiTo+C*lnCC-lnPhiT) / delta_mua;
  else
    clear foo;
    foo = FwdModel.P(idxLambda).PhiTotal;
    lnPhiTo(1:nMeas) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnPhiTo(nMeas+1:nMeas*2) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
    lnCC = -(C'*C) \ (C' * (lnPhiTo - lnPhiM));
    dlnPhiT_dmua = (lnPhiTo+C*lnCC-lnPhiT) / delta_mua;
  end

  % create the Jacobian
  clear J
  J(:,1) = dlnPhiT_dmusp;
  J(:,2) = dlnPhiT_dmua;
  
  s = diag(1./sum(J.*J).^0.5);
  J = J * s;

  % get the difference
  y = lnPhiT - lnPhiM;
  lsq = sum(y.^2);
  
  % calculate the step
  delta_mu = - s * ((J' * J + lambda* eye(size(J,2))) \ (J' * y));
  
  % calculate CCs and CCd
%  CCs(1) = 1;
%  CCs(2:nSrcs,1) = exp(lnCCo(1:nSrcs-1));
%  CCd = exp(lnCCo(nSrcs:nCC));
  
  CCs = zeros(length(pmiModel.Src.Pos),1);
  CCd = zeros(length(pmiModel.Det.Pos),1);

  if pmiModel.ModFreq == 0
    for idx = 1:length(pmiModel.Src.Pos)
      clear foo;
      foo = find( idx==pmiModel.MeasList(idxThisLambda,1) );
      if ~isempty( foo )
	if srcIdx(idx)~=1
	  CCs(idx) = exp(lnCCo(srcIdx(idx)-1));
	else
	  CCs(idx) = 1;
	end
      end
    end
    for idx = 1:length(pmiModel.Det.Pos)
      clear foo;
      foo = find( idx==pmiModel.MeasList(idxThisLambda,2) );
      if ~isempty( foo )
	CCd(idx) = exp(lnCCo(nSrcs+detIdx(idx)-1));
      end
    end
  else
    for idx = 1:length(pmiModel.Src.Pos)
      clear foo;
      foo = find( idx==pmiModel.MeasList(idxThisLambda,1) );
      if ~isempty( foo )
	if srcIdx(idx)~=1
	  CCs(idx) = exp(lnCCo(srcIdx(idx)-1) + sqrt(-1)* ...
			 lnCCo(nCC+srcIdx(idx)-1) );
	else
	  CCs(idx) = 1;
	end
      end
    end
    for idx = 1:length(pmiModel.Det.Pos)
      clear foo;
      foo = find( idx==pmiModel.MeasList(idxThisLambda,2) );
      if ~isempty( foo )
	CCd(idx) = exp(lnCCo(nSrcs+detIdx(idx)-1) + sqrt(-1)* ...
		       lnCCo(nCC+nSrcs+detIdx(idx)-1) );
      end
    end
  end