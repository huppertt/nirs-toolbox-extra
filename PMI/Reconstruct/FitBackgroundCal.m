%function [Fmu, lsqo, CCs, CCd, P]=FitBackgroundCal(Mu,ds,param);
% Mu is the initial guess for the unknowns we are solving for.
% This routine fits for the background optical properties and
% the calibration coefficients.
%
% (Mu_sp, Mu_a) for each wavelength.
%
%  Mu_sp(lambda1) = Mu(1,1);
%  Mu_sp(lambda2) = Mu(1,2);
%
%  Mu_a(lambda1) = Mu(2,1);
%  Mu_a(lambda2) = Mu(2,2);
%
%
%  
%
% ds is the data structure containing the experimental data
%  we are fitting ds.PhiTotalN for mu_sp and mu_a
%
% This routine returns the best fit for the mu_sp and mu_a
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/04/19 19:47:35 $
%
%  $Revision: 1.4 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fmu, lsqo_lambda, CCs, CCd, P]=FitBackgroundCal(Mu,ds,param);
  
%
% Do a Newton-Raphson fit
%
  lsqo_lambda = zeros(length(ds.Inv.Lambda),1);
  for idxLambda = 1:length(ds.Inv.Lambda)
    foo = find(idxLambda==ds.Inv.MeasList(:,4));
    if ~isempty(foo)
      lambda = param.lambda;
      niter = 1;
      
      [delta_mu lsqo CCs(:,idxLambda) CCd(:,idxLambda) P(idxLambda).PhiInc] = ...
	  FitBackgroundCal1(Mu(:,idxLambda), ds.Inv, idxLambda, ...
			    lambda);
      
      if param.verbose
	[Mu(:,idxLambda) delta_mu]
	lsqo
      end
      
      Mu_temp(:,idxLambda) = Mu(:,idxLambda) + delta_mu;
      while( lambda<param.lambdaMax & niter<10000)
	niter = niter + 1;
	
	[delta_mu lsq CCs_temp CCd_temp PhiInc_temp] = ...
	    FitBackgroundCal1(Mu_temp(:,idxLambda), ds.Inv, idxLambda, ...
			      lambda);
	
	
	if param.verbose
	  [Mu_temp(:,idxLambda) delta_mu]
	  lsq
	end
	
	if( lsq<lsqo)
	  Mu(:,idxLambda) = Mu_temp(:,idxLambda);
	  CCs(:,idxLambda) = CCs_temp;
	  CCd(:,idxLambda) = CCd_temp;
	  P(idxLambda).PhiInc = PhiInc_temp;
	  Mu_temp(:,idxLambda) = Mu_temp(:,idxLambda) + delta_mu;
	  lsqo = lsq;
	  lambda = lambda / param.lambdaDec;
	else
	  lambda = lambda * param.lambdaInc;
	end
	
	
	if param.verbose
	  lambda
	end
	
      end
    end
  
    lsqo_lambda(idxLambda,1) = lsqo;
  end
  Fmu = Mu;


return

  
