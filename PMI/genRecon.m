%GenRecon       Generate the selected Reconstruction.
%
%   pmi = genRecon(pmi)
%
%   pmi          The Photon Migration Imaging data structure with result
%                included.
%
%   GenRecon generates the reconstruction specified in the PMI data structure
%   and if requested prints out performance measures of the reconstruction.
%
%   Calls: art, sirt, fatmn, fattsvd, fatmtsvd, tcgls
%
%   Bugs: none known.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%  $Author: dboas $
%
%  $Date: 2001/08/16 23:02:37 $
%
%  $Revision: 1.10 $
%
%  $Log: genRecon.m,v $
%  Revision 1.10  2001/08/16 23:02:37  dboas
%  Fixed a minus sign for reconstructing delta_musp
%
%  Revision 1.9  2001/05/18 18:49:43  dboas
%  Changed reconstruction of delta_musp to be -delta_musp / musp_o rather than
%  delta_D / Do.
%
%  Revision 1.8  2001/04/17 18:39:44  dboas
%  Fixed a matrix dimension problem for handle multi-wavelength scattering
%  perturbation.
%
%  Revision 1.7  2001/02/25 20:34:41  dboas
%  No need to scale delta_D recon as the matrix is already scaled for delta_D / D.
%
%  Revision 1.6  2001/02/21 20:22:33  dboas
%  Scaling of the A matrix fixed for multi-wavelength.
%
%  Revision 1.5  2000/11/30 20:34:09  dboas
%  Fixed ds.Recon.Musp to show image of delta_musp rather than delta_D / D_o.
%
%  Revision 1.4  2000/10/12 18:08:49  dboas
%  Scaling matrix by muao and muspo to make unitless and to put different unknowns
%  on equal footing.
%
%  TSVDnSV changed to TSVD_nSV
%
%  Tikhonov regularization added.
%
%  Revision 1.3  2000/09/01 19:57:07  dboas
%  Added a pmi.Recon.Whiten flag so that the user can specify whether to whiten
%  the matrix (=1) or not (=0).  It does not make sense to whiten the Rytov
%  approximation (at the moment) or Maximum Likelihood.
%
%  We have removed PhiScat from the Model.P structure.  PhiScat is now calculated
%  in genRecon when needed.
%
%  We have removed the whitened variables (Aw, PhiTotalNw, PhiIncw, etc.) from the
%  Model.P structure.  These are now calculated within genRecon when needed.
%
%  Fixed up flags for TSVD a little.  Flags added to pmi.Recon are
%        pmi.Recon.TSVD_FullSVS - sometimes the full SVS does not converge
%  			     and thus setting this to 0 causes only the
%  			     necessary (truncation level) elements of the
%  			     SVS to be calculated.  This approach can be
%  			     slower but has a better chance of converging
%  			     in those rare instances.  Also, setting this
%  			     to FALSE (=0) has the disadvantage of requiring
%  			     additional calls to SVDS if the truncation level
%  			     is increased.  The default is 1.
%        pmi.Recon.TSVD_Lsq - Determines whether to calculate the SVD of
%  			 A'A (the least squares solution =1) or A (=0).
%  			 The default is 0.
%        pmi.Recon.TSVD_CalcSVD - if true (=1) then recalculate the SVD otherwise
%  			     use the previously calculated SVD.  No default,
%  			     this variable must be set by the user.
%
%  We call fattsvd with the options.  fattsvdls has been removed from the package
%  as it is redundant with fattsvd.
%
%  Revision 1.2  2000/07/27 15:00:23  dboas
%  Added P(idxLambda) sub-structure.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.2  2000/02/12 00:12:36  dboas
%  Use pmi.PhiTotalNw for determining nLambda instead of pmi.Inv.Scatw
%
%  Revision 3.1  2000/02/02 19:40:55  rjg
%  Scattered field estimate is now expected to exist.
%
%  Revision 3.0  1999/06/17 19:29:38  rjg
%  Initial Revision for PMI 3.0
%
%  Revision 2.1  1999/02/05 20:43:34  rjg
%  Added help brief comments.
%
%  Revision 2.0  1999/02/05 20:40:09  rjg
%  *** empty log message ***
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmi = genRecon(pmi)

%
% Copy A into local variable and whiten if necessary
% Calculate PhiScat and put in local variable and whiten if necessary
%
  if ~isfield(pmi.Recon,'Whiten')
    pmi.Recon.Whiten = 0;
  end
  
  if pmi.Recon.Whiten & (strcmpi(pmi.Inv.Method.Type,'Rytov') ...
	  | strcmpi(pmi.Recon.ReconAlg,'MaxLike'))
    error('You should not whiten the Rytov Approximation or MaxLike')
  end
  
  if strcmp('uniform', lower(pmi.Inv.CompVol.Type))
    nX = length(pmi.Inv.CompVol.X);
    nY = length(pmi.Inv.CompVol.Y);
    nZ = length(pmi.Inv.CompVol.Z);
    nVoxs = nX * nY *nZ;
  else
    nVoxs = size(pmi.Inv.CompVol.Pos);
  end

  if isfield(pmi.Inv.Method,'ObjVec_mua')
    calc_mua = pmi.Inv.Method.ObjVec_mua;
  else
    calc_mua = 0;
  end
  if isfield(pmi.Inv.Method,'ObjVec_musp')
    calc_musp = pmi.Inv.Method.ObjVec_musp;
  else
    calc_musp = 0;
  end
  
  nLambda = length(pmi.Inv.Lambda);
  for idx = 1:nLambda
    if pmi.Recon.Whiten
      pmi.Noise.P(idx).w = 1 ./ ...
	  (pmi.Noise.P(idx).TotalVar).^0.5;
      P(idx).A = rowscale(pmi.Inv.P(idx).A, ...
					 pmi.Noise.P(idx).w);    
      P(idx).PhiScat = (pmi.Inv.P(idx).PhiTotalN - ...
	  pmi.Inv.P(idx).PhiInc) .* pmi.Noise.P(idx).w;
    elseif strcmp(pmi.Inv.Method.Type,'Born')
      P(idx).A = pmi.Inv.P(idx).A;
      P(idx).PhiScat = pmi.Inv.P(idx).PhiTotalN - ...
	  pmi.Inv.P(idx).PhiInc;
    elseif strcmp(pmi.Inv.Method.Type,'Rytov')
      P(idx).A = pmi.Inv.P(idx).A;
      P(idx).PhiScat = Complex_Log_PhiByPhi(pmi.Inv, ...
			pmi.Inv.P(idx).PhiTotalN, ...
			pmi.Inv.P(idx).PhiInc, idx);
    end      
  
    %%
    %% Scale Matrix:  The mua unknown is normalized by  Mu_ao
    %%                The delta_D unknown is already normalized by Do
    %%
    if calc_mua
      P(idx).A(:,1:nVoxs) = P(idx).A(:,1:nVoxs) * pmi.Inv.Mu_a(idx);
    end
%    if calc_musp
%      if calc_mua
%	P(idx).A(:,nVoxs+1:2*nVoxs) = P(idx).A(:,nVoxs+1:2*nVoxs) * ...
%	    pmi.Inv.Mu_sp(idx);
%      else
%	P(idx).A(:,1:nVoxs) = P(idx).A(:,1:nVoxs) * pmi.Inv.Mu_sp(idx);
%      end
%    end
  
  end

  

%%
%%  Calculate an esimate of the weighted scattered field
%%
%WScatEst = pmi.PhiTotalNw - pmi.Noise.w .* pmi.Inv.PhiInc;
%nLambda = length(pmi.Inv.Lambda);
%for idx=1:nLambda
%  if strcmpi(pmi.Inv.Method.Type, 'Born')
%    pmi.Inv.P(idx).PhiScatw = pmi.Inv.P(idx).PhiTotalNw - pmi.Inv.P(idx).PhiIncw;
%    pmi.Inv.P(idx).PhiScat = pmi.Inv.P(idx).PhiTotalN - pmi.Inv.P(idx).PhiInc;
%  elseif strcmpi(pmi.Inv.Method.Type, 'Rytov')
%    pmi.Inv.P(idx).PhiScatw = Complex_Log_PhiByPhi(pmi.Inv, ...
%				  pmi.Inv.P(idx).PhiTotalNw, ...
%				  pmi.Inv.P(idx).PhiIncw, idx);
%    pmi.Inv.P(idx).PhiScat = Complex_Log_PhiByPhi(pmi.Inv, ...
%				  pmi.Inv.P(idx).PhiTotalN, ...
%				  pmi.Inv.P(idx).PhiInc, idx);
%    pmi.Inv.P(idx).PhiScatw = log(pmi.Inv.P(idx).PhiTotalNw ./ pmi.Inv.P(idx).PhiIncw);
%    pmi.Inv.P(idx).PhiScat = log(pmi.Inv.P(idx).PhiTotalN ./ pmi.Inv.P(idx).PhiInc);
%  end
%end

%%
%%  Reconstruct image of volume
%%
nUnknowns = size(P(1).A, 2);

%%
%%  Compute x estimate using requested algorithm
%%
xEst = zeros(nUnknowns, nLambda);
switch pmi.Recon.ReconAlg

case 'Back Projection'
    for j = 1:nLambda
        %%
        %%  Rescale the backprojection matrix such test the norm of the rows
        %%  of the transpose of A are equal to 1.
        %%
        BPO = zeros(nUnknowns, size(P(j).A, 1));
        for i = 1:nUnknowns
            temp = P(j).A(:,i)';
            BPO(i,:) = temp ./ norm(temp);
        end
        xEst(:,j) = BPO * P(j).PhiScat;
    end

case 'ART'
    for j = 1:nLambda
        xEst(:,j) = art(P(j).A, P(j).PhiScat, ...
            zeros(nUnknowns, 1), pmi.Recon.ART_nIter);
    end

case 'SIRT'
    for j = 1:nLambda
        xEst(:,j) = sirt(P(j).A, P(j).PhiScat, ...
            zeros(nUnknowns, 1), pmi.Recon.SIRT_nIter);
    end
    
case 'Min. Norm'
    for j = 1:nLambda
        xEst(:,j) = fatmn(P(j).A, P(j).PhiScat);
    end
    
 case 'TSVD'
  option.FullSVS = 1;
  option.Lsq = 0;
  if isfield(pmi.Recon,'TSVD_FullSVS')
    option.FullSVS = pmi.Recon.TSVD_FullSVS;
  end
  if isfield(pmi.Recon,'TSVD_Lsq')
    option.Lsq = pmi.Recon.TSVD_Lsq;
  end
  
  for j = 1:nLambda
    if pmi.Recon.TSVD_CalcSVD
      if pmi.Debug
	disp('Computing SVD');
      end
      
      [xEst(:,j) U S V Ainv] = ...
	  fattsvd(P(j).A, P(j).PhiScat, ...
		    pmi.Recon.TSVD_nSV, option);
      pmi.Recon.W(j).Uecon = U;
      pmi.Recon.W(j).Secon = S;
      pmi.Recon.W(j).Vecon = V;
      pmi.Inv.P(j).Ainv = Ainv;

      
      if j == nLambda
	pmi.Recon.TSVD_CalcSVD = 0;
      end
    else
      if pmi.Debug
	disp('Already have SVD');
      end
      [xEst(:,j) foo1 foo2 foo3 Ainv] = fattsvd(P(j).A, P(j).PhiScat, ...
			  pmi.Recon.TSVD_nSV, option, pmi.Recon.W(j).Uecon, ...
			  pmi.Recon.W(j).Secon, pmi.Recon.W(j).Vecon);
      pmi.Inv.P(j).Ainv = Ainv;
    end
    
  end
    
case 'MTSVD'
    LaplOp = lapl3d(length(pmi.Inv.CompVol.X),length(pmi.Inv.CompVol.Y),...
        length(pmi.Inv.CompVol.Z));
    for j = 1:nLambda
        if pmi.Recon.TSVD_CalcSVD
            %%
            %%  If the system has changed we must recompute the full SVD as
            %%  well as the MTSVD solution
            %%
            if pmi.Debug
                fprintf('Computing full MTSVD\n')
            end
            pmi.Recon.xTSVD = zeros(nUnknowns, nLambda);
            pmi.Recon.xNS = zeros(nUnknowns, nLambda);
            pmi.Recon.S = [];
            pmi.Recon.V = [];
            pmi.Recon.U = [];
            [xEst(:,j) pmi.Recon.xTSVD(:,j) pmi.Recon.xNS(:,j) ...
                    pmi.U(:,:,j) pmi.S(:,:,j) pmi.V(:,:,j)] = ...
                fatmtsvd(pmi.Inv.Aw(:,:,j), pmi.Inv.PhiScatw(:,j), ...
                pmi.Recon.MTSVDnSV, LaplOp, pmi.Recon.MTSVDLambda);
            prev_nSV = pmi.Recon.MTSVDnSV;
            if j == nLambda
                pmi.Recon.TSVD_CalcSVD = 0;
            end

        elseif pmi.Recon.MTSVDnSV == prev_nSV
            %%
            %%  If all test has been modified is lambda we can quickly compute
            %%  the result.
            %%
            if pmi.Debug
                fprintf('Only changed lambda\n')
            end
            xEst(:,j) = pmi.Recon.xTSVD(:,j) - ...
                pmi.Recon.MTSVDLambda * pmi.Recon.xNS(:,j);
            
    
        else    
            %%
            %%  If the system is still the same then recompute using the
            %%  existing SVD computation
            %% 
            if pmi.Debug
                fprintf('Recomputing optimum NULL space contribution\n')
            end
            [xEst(:,j) pmi.Recon.xTSVD(:,j) pmi.Recon.xNS(:,j)] = ...
                fatmtsvd(pmi.Inv.Aw(:,:,j), pmi.Inv.PhiScatw(:,j), ...
                pmi.Recon.MTSVDnSV, LaplOp, pmi.MTSVDLambda, ...
                pmi.Recon.U(:,:,j), pmi.Recon.S(:,:,j), pmi.Recon.V(:,:,j));
            prev_nSV = pmi.Recon.MTSVDnSV;
        end
    end

case 'TCG'
    %%
    %%  Truncated normal equation solution
    %%
    for j = 1:nLambda
        xEst(:,j) = tcgls(P(j).A, P(j).PhiScat, ...
            pmi.Recon.TCGnIter);
    end

case 'MaxLike'
    %%
    %%  Maximum likelihood solution
    %%
    for idx = 1:nLambda
      if ~strcmpi(pmi.Inv.Method.Type,'Rytov')
	C = diag(pmi.Noise.P(idx).TotalVar);
      else
	C = diag(pmi.Noise.P(idx).TotalVar./(pmi.Inv.P(idx).PhiInc.^2) ...
		 );
      end
      pmi.Inv.P(idx).Ainv = P(idx).A' * inv(P(idx).A*P(idx).A' + C );
      xEst(:,idx) = pmi.Inv.P(idx).Ainv * P(idx).PhiScat;
    end

case 'Tik'
   %%
   %% Tikhanov Regularization
   %%

   for idx = 1:nLambda
      
      foo = P(idx).A * P(idx).A';
      lambda = pmi.Recon.TIK_alpha * max(max(foo));
      pmi.Inv.P(idx).Ainv = P(idx).A' * inv(foo + lambda*eye(size(foo)));
      xEst(:,idx) = pmi.Inv.P(idx).Ainv * P(idx).PhiScat;

   end
    
    
otherwise
    error(['Unknown Reconstruction technique: ' pmi.Recon.ReconAlg]);
end

nX = length(pmi.Inv.CompVol.X);
nY = length(pmi.Inv.CompVol.Y);
nZ = length(pmi.Inv.CompVol.Z);
nVoxels = nX * nY * nZ;
idx = 1;
if isfield(pmi.Inv.Method,'ObjVec_mua');
  if pmi.Inv.Method.ObjVec_mua == 1
    clear pmi.Recon.Mua;
    pmi.Recon.Mua = xEst(idx:nVoxels,:);
    pmi.Recon.Mua = pmi.Recon.Mua .* (ones(nVoxels,1)*pmi.Inv.Mu_a);
    idx = nVoxels+1;
  end
end
if isfield(pmi.Inv.Method,'ObjVec_musp');
  if pmi.Inv.Method.ObjVec_musp == 1
    clear pmi.Recon.Musp;
    pmi.Recon.Musp = -(ones(nVoxels,1)*pmi.Inv.Mu_sp) .* ...
	(xEst(idx:idx+nVoxels-1,:));
    idx = idx + nVoxels;
  end
end
if isfield(pmi.Inv.Method,'ObjVec_sd');
  if pmi.Inv.Method.ObjVec_sd == 1
    clear pmi.Recon.SD;
    if strcmp(pmi.Inv.Method.Type,'Rytov')
      pmi.Recon.SD = exp( xEst(idx:size(xEst,1),:) );
    else
      pmi.Recon.SD = xEst(idx:size(xEst,1),:);
    end
  end
end


%%
%%  Print out performance measures of the Reconstruction
%%
if pmi.Debug > 0
    pmi.Recon.xError = [];
    pmi.Recon.xResid = [];
    nLambda = size(pmi.Fwd.Mu_a, 2);
    for j = 1:nLambda
        delMu_a = calcDelMuA(pmi.Object, pmi.Inv.CompVol, pmi.Inv.Mu_a, j);
        delMu_a = delMu_a(:);
        pmi.Recon.xError(:,j) = xEst(:,j) - delMu_a;
        fprintf('Lambda #%d\n', j);
        fprintf('  Normalized 2-norm of the error: %f\n', ...
                norm(pmi.Recon.xError(:,j))/norm(delMu_a));
        fprintf('  Mean abs value of the error: %f\n', ...
            mean(abs(pmi.Recon.xError(:,j))));

        fprintf('  Maximum difference %f\n', max(pmi.Recon.xError(:,j)));
        fprintf('  Minimum difference %f\n', min(pmi.Recon.xError(:,j)));

        pmi.Recon.xResid(:,j) = pmi.Inv.Aw(:,:,j) * xEst(:,j) - ...
            pmi.Inv.PhiScatw(:,j);
        fprintf('  2-norm of residual: %f\n', norm(pmi.Recon.xResid(:,j)));
    end
end
