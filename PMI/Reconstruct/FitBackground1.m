%function [delta_mu, lsq, CCs, CCd, PhiInc]=FitBackground1(Mu,pmiModel,idxLambda, lambda)
%
% This is a sub-function to FitBackground and is used to get
% the Newton-Raphson update
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/07/27 15:09:23 $
%
%  $Revision: 1.2 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delta_mu, lsq, CCs, CCd, PhiInc]=FitBackground1(Mu,pmiModel,idxLambda, lambda)

  delta_musp = 0.001;
  delta_mua = 0.00001;
  
  if pmiModel.ModFreq == 0
    lnPhiM = log(pmiModel.P(idxLambda).PhiTotalN);
  else
    nMeas = size(pmiModel.MeasList,1);
    clear foo;
    foo = pmiModel.P(idxLambda).PhiTotalN;
    lnPhiM(1:nMeas,1) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnPhiM(nMeas+1:nMeas*2,1) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
  end
  

  %Evaluate at Mu_o
  pmiModel.Mu_sp(idxLambda) = Mu(1);
  pmiModel.Mu_a(idxLambda) = Mu(2);
%  amp = Mu(3);
  FwdModel=genMeasData(pmiModel);
  PhiInc = FwdModel.P(idxLambda).PhiInc;
  if pmiModel.ModFreq == 0
    lnGT = log(FwdModel.P(idxLambda).PhiTotal);
    lnAmp = mean(lnPhiM - lnGT);
    amp0 = exp(lnAmp); phase0 = 0;
    lnPhiT = lnAmp + lnGT;
  else
    clear foo;
    foo = FwdModel.P(idxLambda).PhiTotal;
    lnGT(1:nMeas) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnGT(nMeas+1:nMeas*2) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
    lnAmp = mean(lnPhiM(1:nMeas) - lnGT(1:nMeas)');
    Phase = mean(lnPhiM(nMeas+1:nMeas*2) - lnGT(nMeas+1:nMeas*2)');
    amp0 = exp(lnAmp); phase0 = Phase;
    lnPhiT(1:nMeas,1) = lnAmp + lnGT(1:nMeas)';
    lnPhiT(nMeas+1:nMeas*2,1) = Phase + lnGT(nMeas+1:nMeas*2)';
  end

  %Evaluate the derivative versus mu_sp
  pmiModel.Mu_sp(idxLambda) = Mu(1)+delta_musp;
  pmiModel.Mu_a(idxLambda) = Mu(2);
  FwdModel=genMeasData(pmiModel);
  if pmiModel.ModFreq == 0
    lnGT = log(FwdModel.P(idxLambda).PhiTotal);
    lnAmp = mean(lnPhiM - lnGT);
    dlnPhiT_dmusp = (lnAmp+lnGT-lnPhiT) / delta_musp;
  else
    clear foo;
    foo = FwdModel.P(idxLambda).PhiTotal;
    lnGT(1:nMeas) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnGT(nMeas+1:nMeas*2) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
    lnAmp = mean(lnPhiM(1:nMeas) - lnGT(1:nMeas)');
    Phase = mean(lnPhiM(nMeas+1:nMeas*2) - lnGT(nMeas+1:nMeas*2)');
    dlnPhiT_dmusp(1:nMeas,1) = (lnAmp + lnGT(1:nMeas)'-lnPhiT(1:nMeas)) /delta_musp;
    dlnPhiT_dmusp(nMeas+1:nMeas*2,1) = (Phase + lnGT(nMeas+1:nMeas*2)'-lnPhiT(nMeas+1:nMeas*2)) /delta_musp;
  end
  
  %Evaluate the derivative versus mu_a
  pmiModel.Mu_sp(idxLambda) = Mu(1);
  pmiModel.Mu_a(idxLambda) = Mu(2) + delta_mua;
  FwdModel=genMeasData(pmiModel);
  if pmiModel.ModFreq == 0
    lnGT = log(FwdModel.P(idxLambda).PhiTotal);
    lnAmp = mean(lnPhiM - lnGT);
    dlnPhiT_dmua = (lnAmp+lnGT-lnPhiT) / delta_mua;
  else
    clear foo;
    foo = FwdModel.P(idxLambda).PhiTotal;
    lnGT(1:nMeas) = log(abs(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2)));
    lnGT(nMeas+1:nMeas*2) = angle(foo(1:nMeas) + (-1)^0.5 * foo(nMeas+1:nMeas*2));
    lnAmp = mean(lnPhiM(1:nMeas) - lnGT(1:nMeas)');
    Phase = mean(lnPhiM(nMeas+1:nMeas*2) - lnGT(nMeas+1:nMeas*2)');
    dlnPhiT_dmua(1:nMeas,1) = (lnAmp + lnGT(1:nMeas)'-lnPhiT(1:nMeas)) /delta_mua;
    dlnPhiT_dmua(nMeas+1:nMeas*2,1) = (Phase + lnGT(nMeas+1:nMeas*2)'-lnPhiT(nMeas+1:nMeas*2)) /delta_mua;
  end

  % create the Jacobian
  clear J;
  J(:,1) = dlnPhiT_dmusp;
  J(:,2) = dlnPhiT_dmua;
  
  s = diag(1./sum(J.*J).^0.5);
  J = J * s;

  % get the difference
  y = lnPhiT - lnPhiM;
  lsq = sum(y.^2);
  
  % calculate the step
  delta_mu = - s* ((J' * J + lambda*eye(size(J,2))) \ (J' * y));

  % calculate CCs and CCd
  nSrcs = size(pmiModel.Src.Pos,1);
  nDets = size(pmiModel.Det.Pos,1);
  
  clear CCS; clear CCd;
  CCs(1:nSrcs,1) = ones(nSrcs,1) * amp0 * exp(j*phase0);
  CCd(1:nDets,1) = ones(nDets,1);

  
  
  
