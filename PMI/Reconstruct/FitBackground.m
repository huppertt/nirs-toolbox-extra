%function [Fmu, lsqo, CCs, CCd, PhiInc]=FitBackground(Mu,ds);
% Mu is the initial guess for the unknowns we are solving for
% (Mu_sp, Mu_a, amp) for each wavelength.
%  Mu_sp(lambda1) = Mu(1,1);
%  Mu_sp(lambda2) = Mu(1,2);
%
%  Mu_a(lambda1) = Mu(2,1);
%  Mu_a(lambda2) = Mu(2,2);
%
%  amp(lambda1) = Mu(3,1);
%  amp(lambda2) = Mu(3,2);
%
% ds is the data structure containing the experimental data
%  we are fitting ds.PhiTotalN for mu_sp and mu_a
%
% This routine returns the best fit for the mu_sp and mu_a and amp
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Fmu, lsqo, CCs, CCd, PhiInc]=FitBackground(Mu,ds);

  
%
% Do a Newton-Raphson fit
%
  for idxLambda = 1:length(ds.Fwd.Lambda)
    
    lambda = 1e-1;
    niter = 1;
    [delta_mu lsqo CCs(:,idxLambda) CCd(:,idxLambda) PhiInc] = FitBackground1(Mu(:,idxLambda), ds.Inv, idxLambda, ...
				     lambda);
    
%    Mu(:,idxLambda)
%    lsqo
    
    Mu_temp(:,idxLambda) = Mu(:,idxLambda) + delta_mu;
    while( lambda<1e5 & niter<100)
      niter = niter + 1;
      [delta_mu lsq CCs(:,idxLambda) CCd(:,idxLambda) PhiInc] = FitBackground1(Mu_temp(:,idxLambda), ds.Inv, ...
				      idxLambda, lambda);
      
%      [Mu_temp(:,idxLambda) delta_mu]
%      lsq

      if( lsq<lsqo)
	Mu(:,idxLambda) = Mu_temp(:,idxLambda);
	Mu_temp(:,idxLambda) = Mu_temp(:,idxLambda) + delta_mu;
	lsqo = lsq;
	lambda = lambda / 1.5;
      else
	lambda = lambda * 2;
      end
      
    end
  end
  
  Fmu = Mu;
  

return

  
