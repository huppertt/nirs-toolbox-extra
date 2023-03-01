% genMeasData    Generate the measured data.
%
%   pmiModel = genMeasData(pmiModel, pmiObject, Debug);
%
%   pmiModel      The PMI Model data structure to updated.
%   pmiObject     The Object data structure.
%
%
%   GenMeasData generates the noiseless measured data using the forward method
%   information present in the PMI imaging data structure.
%
%   Calls: none.
%
%   Bugs: CURRENTLY ONLY HANDLES 1st ORDER FOR BORN CASE.

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/08/24 20:33:33 $
%
%  $Revision: 1.13 $
%
%  $Log: genMeasData.m,v $
%  Revision 1.13  2001/08/24 20:33:33  dboas
%  Small typo for Scattering Perturbation in Rytov.
%
%  Revision 1.12  2001/08/16 23:01:29  dboas
%  Rytov was correct to handle RF data correctly.
%  Using delta_D/D instead of -delta_musp/musp.
%
%  Revision 1.11  2001/07/13 20:14:18  dboas
%  Changed the scattering perturbation for Born and Rytov to be delMu_sp /
%  Mu_sp_o.  Also corrected a minus sign error for the scattering perturbation.
%
%  Revision 1.10  2001/04/13 23:11:11  dboas
%  Modified ExactSphere to properly parse result for DC mixed with RF.
%
%  Revision 1.9  2001/02/07 16:26:56  dboas
%  Added forward model method type 'ExtBornN'.  This does an iterative extended
%  Born approximation.
%  DON'T FORGET THAT 'BornN' AND 'ExtBornN' DON'T NECESSARILY CONVERGE IF THE
%  PERTURBATION IS TOO LARGE.
%
%  Added Object Type 'Image'.
%
%  Revision 1.8  2001/02/02 17:51:18  dboas
%  dded Extended Born method.
%  The exact sphere calculation had a major error.  Basically, the calculation of
%  k inside the sphere
%  used D from outside the sphere.  This was fixed.
%
%  Revision 1.7  2000/11/07 16:36:08  dboas
%  Added Exact Sphere calculation
%
%  Revision 1.6  2000/09/01 19:27:20  dboas
%  Having some trouble with re-initializing the P sub-structure.  This is probably
%  still not the desired fix.
%
%  Added Model.Method.Type 'FullBorn' and 'BornN'.  For 'BornN' it is necessary to
%  specify 'Model.Method.Born_Order'
%
%  Revision 1.5  2000/08/01 13:23:13  tgaudett
%  More Fixes for Rytov
%
%  Revision 1.4  2000/08/01 12:07:55  tgaudett
%  Fixed Rytov Approximation Code for P stucture.
%
%  Revision 1.3  2000/07/27 14:55:11  dboas
%  Added the P(idxLambda) sub-structure to handle data of different wavelengths.
%  The optode amplitudes must now be 2 dimensional.
%
%  Revision 1.2  2000/06/27 14:16:40  dboas
%  Added scattering perturbation to the First Born approximation
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.6  2000/04/13 20:16:41  dboas
%  pmi.PhiTotal assignment added to Helmholtz Homogeneous calculation
%
%  Revision 3.5  2000/01/10 00:14:14  dboas
%  Storing the source and detector lists for use by other functions
%
%  Revision 3.4  1999/12/03 13:52:23  dboas
%  Calculate the incident to purturbation ratio only if
%  the PhiScat field has been calculated.
%  This is needed because this routine
%  can now just calculate the incident
%  fluence (previous revision).
%
%  Revision 3.3  1999/11/18 17:37:37  rjg
%  Correct remaining calls to FullMeasList
%
%  Revision 3.2  1999/11/16 22:37:46  rjg
%  Added 'Helmholtz Homogenous' section from David.
%
%  Revision 3.1  1999/11/15 19:13:39  rjg
%  Commented out extended born section since that was moved to the non-core
%  development directory will uncomment when it is working.
%
%  Fixed calculation of effective source position for FDFD method.
%
%  Added warning about spherical harmonic not being fully functional.
%
%  Spherical harmonic code updated to handle new Object format.
%
%  Removed checking for nLambda since Lambda is avaiable in the new structure.
%
%  Revision 3.0  1999/06/17 19:29:38  rjg
%  Initial Revision for PMI 3.0
%
%  Revision 2.2  1999/02/05 20:39:28  rjg
%  Correctly handles extracting number of wavelengths from mu_a
%  Handles both spheres and blocks (multiple)
%  Calculates the total variance of the noise for appropriate weighting
%  Added TotalVar field to structure
%  No need for update flag for noise wieghting, removed.
%  Added code to added vector norm weighted noise (Eric usual manner
%  for specifying SNR).  Does not correctly handle multiple wavelengths yet.
%
%  Revision 2.1  1998/09/09 15:12:21  rjg
%  Corrected scaling for detector noise to be relative to the unpurturbed
%  scattered field.  Not significant unless two or more noise types were used.
%
%  Revision 2.0  1998/08/20 18:58:22  rjg
%  *** empty log message ***
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmiModel = genMeasData(pmiModel, pmiObject, Debug);

%  if ~exist('Debug')  %doesn't work for WINDOWS
    Debug = 0;
%  end
  
%%
%%  Generate measured data
%%

switch pmiModel.Method.Type

case 'Matlab Variable' 
    %%
    %%  Matlab variable data should be a matrix with 2 columns, the first ...
    %%  column should contain the total fluence measurement and the second ...
    %%  column should contain the incident fluence level.
    %%
    if isfield(pmiModel,'PhiTotal')
      pmiModel = rmfield(pmiModel,'PhiTotal');
    end
    if isfield(pmiModel,'PhiInc')
      pmiModel = rmfield(pmiModel,'PhiInc');
    end
    if isfield(pmiModel,'PhiScat')
      pmiModel = rmfield(pmiModel,'PhiScat');
    end
    
    Fluence = eval(pmiModel.Method.MatlabVarName);
    pmi.PhiTotal = Fluence(:,1);
    if pmiModel.ModFreq > 0
        pmiModel.PhiInc = Fluence(:,2);
        pmiModel.PhiScat = pmi.PhiTotal - pmiModel.PhiInc;
    end
    
case 'Helmholtz Homogeneous'

 nLambda = length(pmiModel.Lambda);
 nFreq = length(pmiModel.ModFreq);

 %%
 %%  Loop over wavelength
 %%
 for idxLambda = 1:nLambda

   if isfield(pmiModel,'P')
     if length(pmiModel.P)>=idxLambda
       if isfield(pmiModel.P(idxLambda),'PhiTotal')
	 pmiModel.P(idxLambda).PhiTotal = [];
       end
       if isfield(pmiModel.P(idxLambda),'PhiInc')
	 pmiModel.P(idxLambda).PhiInc = [];
       end
       if isfield(pmiModel.P(idxLambda),'PhiScat')
	 pmiModel.P(idxLambda).PhiScat = [];
       end
     end
   end
 
   %%
   %%  Loop over frequency
   %%
   idxRow = 1;
   for idxFreq = 1:nFreq
     idxThisFreq = find((idxFreq == pmiModel.MeasList(:,3)) & ...
			(idxLambda == pmiModel.MeasList(:,4)));
     MeasList = pmiModel.MeasList(idxThisFreq,:);
     
     if ~isempty(MeasList)
       PhiInc = DPDWHelmholtz( pmiModel, MeasList, Debug);
       
       %%
       %%  Scale the amplitude of the incident field if 
       %%  necessary.
       %%
       clear SDWeight;
       for idx = 1:size(MeasList,1)
	 SDWeight(idx,1) = ...
	     pmiModel.Src.Amplitude(MeasList(idx,1),idxLambda) ...
	     * pmiModel.Det.Amplitude(MeasList(idx,2),idxLambda);
       end
       PhiInc = SDWeight .* PhiInc;
       
       %%
       %%  Copy the fluence into the correct block
       %%
       nMeas = size(MeasList, 1);
       pmiModel.P(idxLambda).PhiInc(idxRow:idxRow+nMeas-1,1) = real(PhiInc);
       idxRow = idxRow + nMeas;            
       if pmiModel.ModFreq(idxFreq) > 0
	 pmiModel.P(idxLambda).PhiInc(idxRow:idxRow+nMeas-1,1) = imag(PhiInc);
	 idxRow = idxRow + nMeas;
       end
     end
   end
   
   pmiModel.P(idxLambda).PhiTotal = pmiModel.P(idxLambda).PhiInc;
   
 end

   

case 'Born'

 %%
 %%  Check to make sure the number of mu_a parameters is the same as the
 %%  number of delMu_a parameters, this is also the number of
 %%  wavelengths.
 %%
 nLambda = size(pmiModel.Mu_a, 2);
 
 %%
 %%  Generate the volume data using the 1st born method
 %%
 if isfield(pmiModel.Method,'ObjVec_mua')
   calc_mua = pmiModel.Method.ObjVec_mua;
 else
   calc_mua = 0;
 end
 if isfield(pmiModel.Method,'ObjVec_musp')
   calc_musp = pmiModel.Method.ObjVec_musp;
 else
   calc_musp = 0;
 end

 if isfield(pmiModel.P(1),'PhiTotal')
   pmiModel.P(1).PhiTotal = [];
 end
 if isfield(pmiModel.P(1),'PhiScat')
   pmiModel.P(1).PhiScat = [];
 end
 for iLambda = 1:nLambda
   
   nMeas = size(pmiModel.P(iLambda).A, 1);
   pmiModel.P(iLambda).PhiScat = zeros(nMeas, 1);

    if calc_mua
      if strcmpi(pmiObject{1}.Type,'Image') 
	delMu_a = pmiObject{1}.Mu_a(:,iLambda) - pmiModel.Mu_a(iLambda);
      else
	delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
			     iLambda);
      end
    end
    if calc_musp
      if strcmpi(pmiObject{1}.Type,'Image') 
	delMu_sp = pmiObject{1}.Mu_sp(:,iLambda) - pmiModel.Mu_sp(iLambda);
      else
	delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, ...
			       pmiModel.Mu_sp, iLambda);
      end
      Mu_sp_o = pmiModel.Mu_sp(iLambda);
    end

%   if calc_mua
%     delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
%			  iLambda);
%   end
%   if calc_musp
%     delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, pmiModel.Mu_sp, ...
%			    iLambda);
%     Mu_sp_o = pmiModel.Mu_sp(iLambda);
%   end

   if calc_mua & calc_musp
     pmiModel.P(iLambda).PhiScat(:) = pmiModel.P(iLambda).A(:,:) * ...
	 [delMu_a(:); Mu_sp_o ./(Mu_sp_o+delMu_sp(:)) - 1 ];
   elseif calc_mua
     pmiModel.P(iLambda).PhiScat(:) = pmiModel.P(iLambda).A * delMu_a(:);
   elseif calc_musp
     pmiModel.P(iLambda).PhiScat(:) = pmiModel.P(iLambda).A(:,:) * ...
	 (Mu_sp_o ./(Mu_sp_o+delMu_sp(:)) - 1);
   end

   pmiModel.P(iLambda).PhiTotal = pmiModel.P(iLambda).PhiInc + ...
       pmiModel.P(iLambda).PhiScat;
        
 end


case 'Rytov'
 if isfield(pmiModel.Method,'ObjVec_mua')
   calc_mua = pmiModel.Method.ObjVec_mua;
 else
   calc_mua = 0;
 end
 if isfield(pmiModel.Method,'ObjVec_musp')
   calc_musp = pmiModel.Method.ObjVec_musp;
 else
   calc_musp = 0;
 end


 %%
 %%  Check to make sure the number of mu_a parameters is the same as the
 %%  number of delMu_a parameters, this is also the number of
 %%  wavelengths.
 %%
 nLambda = size(pmiModel.Mu_a, 2);
 
 %%
 %%  Generate the volume data using the 1st born method
 %%
 if isfield(pmiModel,'P')
	 if isfield(pmiModel.P(1),'PhiTotal')
     		 pmiModel.P(1).PhiTotal = [];
 	end
 	if isfield(pmiModel.P(1),'PhiScat')
      		pmiModel.P(1).PhiScat = [];
	 end
 end;
 for iLambda = 1:nLambda
   
    nMeas = size(pmiModel.P(iLambda).A, 1);
    pmiModel.P(iLambda).PhiScat = zeros(nMeas,1);
    if calc_mua
      if strcmpi(pmiObject{1}.Type,'Image') 
	delMu_a = pmiObject{1}.Mu_a(:,iLambda) - pmiModel.Mu_a(iLambda);
      else
	delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
			     iLambda);
      end
    end
    if calc_musp
      if strcmpi(pmiObject{1}.Type,'Image') 
	delMu_sp = pmiObject{1}.Mu_sp(:,iLambda) - pmiModel.Mu_sp(iLambda);
      else
	delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, ...
			       pmiModel.Mu_sp, iLambda);
      end
      Mu_sp_o = pmiModel.Mu_sp(iLambda);
    end

   if calc_mua & calc_musp
     pmiModel.P(iLambda).PhiScat(:) = pmiModel.P(iLambda).A(:,:) * ...
	 [delMu_a(:); Mu_sp_o ./(Mu_sp_o+delMu_sp(:)) - 1 ];
   elseif calc_mua
     pmiModel.P(iLambda).PhiScat(:) = pmiModel.P(iLambda).A * delMu_a(:);
   elseif calc_musp
     pmiModel.P(iLambda).PhiScat(:) = pmiModel.P(iLambda).A(:,:) * ...
	 ( Mu_sp_o ./(Mu_sp_o+delMu_sp(:)) - 1);
   end

    NN=size(pmiModel.P(iLambda).PhiScat,1);
    complexPhiInc=pmiModel.P(iLambda).PhiInc(1:NN/2,:)+sqrt(-1)*pmiModel.P(iLambda).PhiInc((NN/2+1):NN,:);
    complexPhiScat=pmiModel.P(iLambda).PhiScat(1:NN/2,:)+sqrt(-1)*pmiModel.P(iLambda).PhiScat((NN/2+1):NN,:);
    
    complexPhiTotal = complexPhiInc .* exp(complexPhiScat);
    pmiModel.P(iLambda).PhiTotal(1:NN/2,:)=real(complexPhiTotal);
    pmiModel.P(iLambda).PhiTotal((NN/2+1):NN,:)=imag(complexPhiTotal);
      
 end


case {'BornN','ExtBornN'}

 nLambda = size(pmiModel.Mu_a, 2);
 
 if isfield(pmiModel.Method,'ObjVec_mua')
   calc_mua = pmiModel.Method.ObjVec_mua;
 else
   calc_mua = 0;
 end
 if isfield(pmiModel.Method,'ObjVec_musp')
   calc_musp = pmiModel.Method.ObjVec_musp;
 else
   calc_musp = 0;
 end

 if isfield(pmiModel.P(1),'PhiTotal')
   pmiModel.P(1).PhiTotal = [];
 end
 if isfield(pmiModel.P(1),'PhiScat')
   pmiModel.P(1).PhiScat = [];
 end
 for iLambda = 1:nLambda

   if calc_mua
     if strcmpi(pmiObject{1}.Type,'Image') 
       delMu_a = pmiObject{1}.Mu_a(:,iLambda) - pmiModel.Mu_a(iLambda);
     else
       delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
			    iLambda);
     end
     nVoxs = length(delMu_a(:));
   end
   if calc_musp
     if strcmpi(pmiObject{1}.Type,'Image') 
       delMu_sp = pmiObject{1}.Mu_sp(:,iLambda) - pmiModel.Mu_sp(iLambda);
     else
       delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, ...
			      pmiModel.Mu_sp, iLambda);
     end
     Mu_sp_o = pmiModel.Mu_sp(iLambda);
     nVoxs = length(delMu_sp(:));
   end
   
%   if calc_mua
%     delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
%			  iLambda);
%     nVoxs = length(delMu_a(:));
%   end
%   if calc_musp
%     delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, pmiModel.Mu_sp, ...
%			    iLambda);
%     Mu_sp_o = pmiModel.Mu_sp(iLambda);
%     nVoxs = length(delMu_sp(:));
%   end

   Pert = pmiModel.P(iLambda).PhiInc_VoxDet;
   PhiTot = Pert;
   P = pmiModel.P(iLambda).A;
   
   if calc_mua & calc_musp
     P(:,1:2*nVoxs) = rowscale(P(:,1:2*nVoxs)',[delMu_a(:); Mu_sp_o ./(Mu_sp_o+delMu_sp(:)) - 1 ])';
     P(:,2*nVoxs+1:3*nVoxs) = rowscale(P(:,2*nVoxs+1:3*nVoxs)', Mu_sp_o ./(Mu_sp_o+delMu_sp(:))-1 )';
     P(:,3*nVoxs+1:4*nVoxs) = rowscale(P(:,3*nVoxs+1:4*nVoxs)', Mu_sp_o ./(Mu_sp_o+delMu_sp(:))-1 )';
   elseif calc_mua
     P(:,1:nVoxs) = rowscale(P(:,1:nVoxs)', delMu_a(:) )';
   elseif calc_musp
     P(:,1:nVoxs) = rowscale(P(:,1:nVoxs)', Mu_sp_o ./(Mu_sp_o+delMu_sp(:))-1 )';
     P(:,nVoxs+1:2*nVoxs) = rowscale(P(:,nVoxs+1:2*nVoxs)', Mu_sp_o ./(Mu_sp_o+delMu_sp(:))-1 )';
     P(:,2*nVoxs+1:3*nVoxs) = rowscale(P(:,2*nVoxs+1:3*nVoxs)', Mu_sp_o ./(Mu_sp_o+delMu_sp(:))-1 )';
   end

   for idx=1:pmiModel.Method.Born_Order
     if calc_musp
       foo = Pert_Grad(Pert, pmiModel, calc_mua);
       Pert = P * foo;
     else
       Pert = P * Pert;
     end
     PhiTot = PhiTot + Pert;
   end
   
   %CONVERT TO MEASLIST ORDER
   Pert = pmiModel.P(iLambda).PhiInc_VoxDet;
   foo = find(pmiModel.MeasList(:,4) == iLambda);
   MeasList = pmiModel.MeasList(foo,:);
   for idx=1:size(MeasList,1)
     PhiTotal(idx,1) = PhiTot(nVoxs+MeasList(idx,2),MeasList(idx, ...
						  1));
     PhiInc(idx,1) = Pert(nVoxs+MeasList(idx,2),MeasList(idx,1));
   end

   %SEPARATE REAL AND IMAGINARY
   if pmiModel.ModFreq == 0
     pmiModel.P(iLambda).PhiTotal = PhiTotal;
     pmiModel.P(iLambda).PhiInc = PhiInc;
   else
     pmiModel.P(iLambda).PhiTotal = [real(PhiTotal); imag(PhiTotal)];
     pmiModel.P(iLambda).PhiInc = [real(PhiInc); imag(PhiInc)];
   end
        
 end


case {'FullBorn','ExtBorn'}

 nLambda = size(pmiModel.Mu_a, 2);
 
 if isfield(pmiModel.Method,'ObjVec_mua')
   calc_mua = pmiModel.Method.ObjVec_mua;
 else
   calc_mua = 0;
 end
 if isfield(pmiModel.Method,'ObjVec_musp')
   calc_musp = pmiModel.Method.ObjVec_musp;
 else
   calc_musp = 0;
 end

 if isfield(pmiModel.P(1),'PhiTotal')
   pmiModel.P(1).PhiTotal = [];
 end
 if isfield(pmiModel.P(1),'PhiScat')
   pmiModel.P(1).PhiScat = [];
 end

 for iLambda = 1:nLambda
   
%   nMeas = size(pmiModel.P(iLambda).A, 1);
%   pmiModel.P(iLambda).PhiScat = zeros(nMeas, 1);

   if calc_mua
     if strcmpi(pmiObject{1}.Type,'Image') 
       delMu_a = pmiObject{1}.Mu_a(:,iLambda) - pmiModel.Mu_a(iLambda);
     else
       delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
			    iLambda);
     end
     nVoxs = length(delMu_a(:));
   end
   if calc_musp
     if strcmpi(pmiObject{1}.Type,'Image') 
       delMu_sp = pmiObject{1}.Mu_sp(:,iLambda) - pmiModel.Mu_sp(iLambda);
     else
       delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, ...
			      pmiModel.Mu_sp, iLambda);
     end
     Mu_sp_o = pmiModel.Mu_sp(iLambda);
     nVoxs = length(delMu_sp(:));
   end
   
%   if calc_mua
%     delMu_a = calcDelMuA(pmiObject, pmiModel.CompVol, pmiModel.Mu_a, ...
%			  iLambda);
%     nVoxs = length(delMu_a(:));
%   end
%   if calc_musp
%     delMu_sp = calcDelMuSp(pmiObject, pmiModel.CompVol, pmiModel.Mu_sp, ...
%			    iLambda);
%     Mu_sp_o = pmiModel.Mu_sp(iLambda);
%     nVoxs = length(delMu_sp(:));
%   end

   Pert = pmiModel.P(iLambda).PhiInc_VoxDet;
   PhiTot = Pert;

   nDets = size(Pert,1) - nVoxs;
   
   if calc_mua & calc_musp
     P = pmiModel.P(iLambda).A( 1:size(pmiModel.P(iLambda).A,1), 1:2*nVoxs+nDets);
     P(:,1:2*nVoxs) = rowscale(P(:,1:2*nVoxs)',[delMu_a(:); Mu_sp_o ./(Mu_sp_o+delMu_sp(:)) - 1 ])';

   elseif calc_mua
     P = pmiModel.P(iLambda).A( 1:size(pmiModel.P(iLambda).A,1), 1:nVoxs+nDets);
     P(:,1:nVoxs) = rowscale(P(:,1:nVoxs)', delMu_a(:) )';

   elseif calc_musp
     P = pmiModel.P(iLambda).A( 1:size(pmiModel.P(iLambda).A,1), 1:nVoxs+nDets);
     P(:,1:nVoxs) = rowscale(P(:,1:nVoxs)', Mu_sp_o ./(Mu_sp_o+delMu_sp(:))-1 )';
   end

   PhiTot = (eye(size(P))-P) \ Pert;
   
   %CONVERT TO MEASLIST ORDER
   Pert = pmiModel.P(iLambda).PhiInc_VoxDet;
   foo = find(pmiModel.MeasList(:,4) == iLambda);
   MeasList = pmiModel.MeasList(foo,:);
   for idx=1:size(MeasList,1)
     PhiTotal(idx,1) = PhiTot(nVoxs+MeasList(idx,2),MeasList(idx,1));
     PhiInc(idx,1) = Pert(nVoxs+MeasList(idx,2),MeasList(idx,1));
   end

   %SEPARATE REAL AND IMAGINARY
   if pmiModel.ModFreq == 0
     pmiModel.P(iLambda).PhiTotal = PhiTotal;
     pmiModel.P(iLambda).PhiInc = PhiInc;
   else
     pmiModel.P(iLambda).PhiTotal = [real(PhiTotal); imag(PhiTotal)];
     pmiModel.P(iLambda).PhiInc = [real(PhiInc); imag(PhiInc)];
   end
        
 end


case 'ExtBorn'
    warning('Extended Born approximation not yet functional');
%    %%
%    %%  Check to make sure the number of Mu_a parameters is the same as the
%    %%  number of delMu_a parameters, this is also the number of
%    %%  wavelengths.
%    %%
%    nLambda = size(pmiModel.Mu_a, 2);
%    nDMu_a = size(pmiModel.Mu_a, 2);
%    if nLambda ~= nDMu_a
%        error('Incorrect number of Mu_a parameters, must match Mu_a');
%    end
%
%    %%
%    %%  Source modulation frequency and source & detector positions
%    %%
%    [pSrc nSrc]= getOptodePos(pmiModel.Src);
%    [pDet nDet]= getOptodePos(pmiModel.Det);
%    nSDPair = nDet * nSrc;
%    nFreq = length(pmiModel.ModFreq);
%
%    nMeas = nSDPair * (nFreq + sum(pmiModel.ModFreq > 0));
%
%    pmiModel.PhiScat = zeros(nMeas, nLambda);
%    pmiModel.PhiInc = zeros(nMeas, nLambda);
%    pmi.PhiTotal = zeros(nMeas, nLambda);
%    if nFreq > 1
%        warning(['Multiple frequencies are not currently implemented for ' ...
%                'Extended born']);
%    end
%
%    %%
%    %%  Create the MeasList matrix
%    %%
%    if ~isfield(pmiModel, 'MeasList')
%        pmiModel.MeasList = FullMeasList(nSrc, nDet, nLambda, nLambda);
%    end
%    
%    %%
%    %%  Calculate the fluence at the detector using the extended born method
%    %%  looping over each wavelength
%    %%
%    for iLambda = 1:nLambda
%        delMu_a = calcDelMuA(pmi.Object, pmiModel.CompVol, pmiModel.Mu_a, iLambda)
%
%        %%
%        %%  Select the appropriate method for the boundary condition
%        %%
%        switch lower(pmiModel.Boundary.Geometry)
%         case { 'semi-infinite', 'semi', 'extrapolated'}
%          if pmi.Debug
%              fprintf(['Executing semi-infinite boundary computation\n']);
%          end
%          %%
%          %%  Move the effective source position 1 mean free path into the medium
%          %%
%          EffpSrc = [pSrc(:,1) pSrc(:,2) pSrc(:,3)-(1/pmiModel.Mu_sp(iLambda))]
%
%          [PhiTotal PhiInc] = DPDWEBornZB(pmiModel.CompVol, pmiModel.Mu_sp(iLambda), ...
%                                          pmiModel.Mu_a(iLambda), delMu_a, ...
%                                          pmiModel.v(iLambda), ...
%                                          pmiModel.idxRefr(iLambda), ...
%                                          pmiModel.ModFreq * 1E6, EffpSrc, pDet, ...
%                                          pmi.Debug);
%
%         case {'infinite', 'inf'}
%          if pmi.Debug
%              fprintf(['Executing infinite boundary computation\n']);
%          end
%          [PhiTotal PhiInc] = DPDWEBornNB(pmiModel.CompVol, pmiModel.Mu_sp(iLambda), ...
%                                          pmiModel.Mu_a(iLambda), delMu_a, ...
%                                          pmiModel.v(iLambda), ...
%                                          pmiModel.idxRefr(iLambda), ...
%                                          pmiModel.ModFreq * 1E6, pSrc, pDet, ...
%                                          pmi.Debug);
%        
%        otherwise
%         error(['Unknown boundary condition: ' pmiFwd.Boundary.Geometry]);
%        end
%
%        if pmiModel.ModFreq > 0
%            pmi.PhiTotal(:, iLambda) = [real(PhiTotal); imag(PhiTotal)];
%            pmiModel.PhiInc(:, iLambda) = [real(PhiInc); imag(PhiInc)];            
%        else
%            pmi.PhiTotal(:, iLambda) = real(PhiTotal);
%            pmiModel.PhiInc(:, iLambda) = real(PhiInc);            
%        end
%    end
%    pmiModel.PhiScat = pmi.PhiTotal - pmiModel.PhiInc;
%
    
case 'FDFD'
    %%
    %%  Adujust the computational volume for the boundary.
    %%  Assuming a negative Z medium
    %%
    if any(strcmp(lower(pmiModel.Boundary.Geometry), {'semi-infinite', 'semi', ...
                    'extrapolated'}))
        %%
        %%  Get the extrapolated boundary distance
        %%
        zBnd = calcExtBnd(pmiModel.idxRefr, pmiModel.Mu_sp);
        if Debug
            fprintf('Extrapolated boundary = %f cm\n', zBnd);
        end
        
        CVShift = - max(pmiModel.CompVol.Z) + zBnd;
        pmiModel.CompVol.Z = pmiModel.CompVol.Z + CVShift;
        if Debug
            fprintf('Z Domain computational volume shifted %f cm\n', CVShift);
        end
    end
    
    %%
    %%  Source modulation frequency and source & detector positions
    %%
    [pSrc nSrc]= getOptodePos(pmiModel.Src);
    [pDet nDet]= getOptodePos(pmiModel.Det);
    nSDPair = nSrc * nDet;
    nFreq = length(pmiModel.ModFreq);
    nMeas = nSDPair * (nFreq + sum(pmiModel.ModFreq > 0));

    nLambda = length(pmiModel.Lambda);
    pmiModel.PhiScat = zeros(nMeas, nLambda);
    pmiModel.PhiInc = zeros(nMeas, nLambda);
    pmi.PhiTotal = zeros(nMeas, nLambda);
    if nFreq > 1
        warning(['Multiple frequencies are not currently implemented for ' ...
                'FDFD']);
    end

    %%
    %%  Create the MeasList matrix
    %%
    if ~isfield(pmiModel, 'MeasList')
        pmiModel.MeasList = FullMeasList(nSrc, nDet, nLambda, nLambda);
    end
    
    %%
    %%  Calculate the fluence at the detector using the extended born method
    %%  looping over each wavelength
    %%
    for iLambda = 1:nLambda
        delMu_a = calcDelMuA(pmi.Object, pmiModel.CompVol, pmiModel.Mu_a, iLambda);

        %%
        %%  Calculate the effiective source position
        %%
        EffpSrc = [pSrc(:,1) pSrc(:,2) pSrc(:,3)-(1/pmiModel.Mu_sp(iLambda))]        
        %%
        %%  Select the appropriate method for the boundary condition
        %%
        switch lower(pmiModel.Boundary.Geometry)
         case { 'semi-infinite', 'semi', 'extrapolated'}
          if Debug
              fprintf(['Executing extrapolated zero boundary' ...
                       ' computation\n']);
          end
          [PhiScat PhiInc] = DPDWFDJacZB(pmiModel.CompVol,  pmiModel.Mu_sp(iLambda), ...
                                         pmiModel.Mu_a(iLambda), delMu_a, ...
                                         pmiModel.v(iLambda),  pmiModel.ModFreq*1E6, ...
                                         EffpSrc, pmiModel.Src.Amplitude,  ...
                                         pDet, pmiModel.Det.Amplitude, Debug);
          

         case {'infinite', 'inf'}
          if Debug
              fprintf(['Executing infinite medium boundary' ...
                       ' computation\n']);
          end
          [PhiScat PhiInc] = DPDWFDJacNB(pmiModel.CompVol,  pmiModel.Mu_sp(iLambda), ...
                                         pmiModel.Mu_a(iLambda), delMu_a, ...
                                         pmiModel.v(iLambda),  pmiModel.ModFreq*1E6, ...
                                         EffpSrc, pmiModel.Src.Amplitude,  ...
                                         pDet, pmiModel.Det.Amplitude, Debug);

         otherwise
          error(['Unknown boundary condition: ' pmiModel.Boundary.Geometry]);
        end

        if pmiModel.ModFreq > 0
            pmiModel.PhiScat(:, iLambda) = [real(PhiScat); imag(PhiScat)];
            pmiModel.PhiInc(:, iLambda) = [real(PhiInc); imag(PhiInc)];
        else
            pmiModel.PhiScat(:, iLambda) = real(PhiScat);
            pmiModel.PhiInc(:, iLambda) = real(PhiInc);            
        end
    end
    pmi.PhiTotal = pmiModel.PhiScat + pmiModel.PhiInc;

    
    
 case 'Spherical'
  
  if ~strcmpi(pmiModel.Boundary.Geometry,'infinite')
    error(['Spherical Exact Forward Solution only works in infinite' ...
	   ' media']);
  end
  
  %%
  %%  Check the object(s) to make sure that they are spheres
  %%
  nObjects = length(pmiObject);
  if nObjects >1
    fprintf(['WARNING: Spherical harmonic solution for' ...
	     ' multiple spheres is\n']);
    fprintf(['only a simple sum approximation of the' ...
	     ' individual scattered fields\n']);
  end
  for iObj = 1:nObjects
    if ~strcmpi(pmiObject{iObj}.Type, 'sphere')
      error('Spherical Harmonic forward solution is only applicable to spherical objects');
    end
  end
  
  %%
  %%
  %%
  nFreq = size(pmiModel.ModFreq,2);
  nLambda = size(pmiModel.Lambda,2);

  %%
  %% Calculate k and D for each wavelength and frequency and object
  %%
  for idxLambda = 1:nLambda
    muspo = pmiModel.Mu_sp(idxLambda);
    muao = pmiModel.Mu_a(idxLambda);
    vo = pmiModel.v(idxLambda);
    v_Do(idxLambda) = vo / (3*muspo);
    for idxFreq = 1:nFreq
      nMeasLF(idxLambda,idxFreq) = length(find(pmiModel.MeasList(:, ...
						  4)==idxLambda & ...
					       pmiModel.MeasList(:, ...
						  3)==idxFreq));
      
      w = pmiModel.ModFreq(idxFreq) * 1e6 * 2 * pi;
      v_ko(idxLambda,idxFreq) = sqrt((-vo*muao+j*w)/v_Do(idxLambda));
      for idxObj = 1:nObjects
	muspi = pmiObject{idxObj}.Mu_sp(idxLambda);
	muai = pmiObject{idxObj}.Mu_a(idxLambda);
	v_Di(idxLambda,idxObj) = vo / (3*muspi);
	v_ki(idxLambda,idxFreq,idxObj) = sqrt((-vo*muai+j*w)/v_Di(idxLambda,idxObj));
      end
    end
  end

  
  %%
  %%  Compute the infinite medium solution for each measurement
  %%
  for idxLambda = 1:nLambda
    
    List = find(pmiModel.MeasList(:,4)==idxLambda);
    nMeas = length(List);

    for idxMeas = 1:nMeas
      idxSrc = pmiModel.MeasList(List(idxMeas),1);
      idxDet = pmiModel.MeasList(List(idxMeas),2);
      rs = pmiModel.Src.Pos(idxSrc,:);
      rd = pmiModel.Det.Pos(idxDet,:);
      idxFreq = pmiModel.MeasList(List(idxMeas),3);
    
      phi_sc = 0;
      for idxObj = 1:nObjects
      
	ro = pmiObject{idxObj}.Pos;
	a = pmiObject{idxObj}.Radius;
	[phi_inc foo] = ExactSphere( rs, ro, rd, v_ko(idxLambda,idxFreq), ...
				     v_Do(idxLambda), v_ki(idxLambda,...
						  idxFreq,idxObj), ... 
				     v_Di(idxLambda,idxObj), a, vo);
	phi_sc = phi_sc + foo;
	
      end

      SD = pmiModel.Src.Amplitude(idxSrc,idxLambda) * ...
	   pmiModel.Det.Amplitude(idxDet,idxLambda);
      PhiTotal(idxMeas,1) = SD*(phi_inc+phi_sc);
      PhiInc(idxMeas,1) = SD*phi_inc;
    end    

    %SEPARATE REAL AND IMAGINARY
    idxMeas = 1;
    for idxFreq = 1:length(pmiModel.ModFreq)
      nMeas = nMeasLF(idxLambda,idxFreq);
      if pmiModel.ModFreq(idxFreq) == 0
	pmiModel.P(idxLambda).PhiTotal(idxMeas:idxMeas+nMeas-1) = ...
	    real(PhiTotal(idxMeas:idxMeas+nMeas-1));
	pmiModel.P(idxLambda).PhiInc(idxMeas:idxMeas+nMeas-1) = ...
	    real(PhiInc(idxMeas:idxMeas+nMeas-1));
	idxMeas = idxMeas + nMeas;
      else
	pmiModel.P(idxLambda).PhiTotal(idxMeas:idxMeas+2*nMeas-1) = ...
	    [real(PhiTotal(idxMeas:idxMeas+nMeas-1));
	     imag(PhiTotal(idxMeas:idxMeas+nMeas-1)) ];
	pmiModel.P(idxLambda).PhiInc(idxMeas:idxMeas+2*nMeas-1) = ...
	    [real(PhiInc(idxMeas:idxMeas+nMeas-1));
	     imag(PhiInc(idxMeas:idxMeas+nMeas-1)) ];
	idxMeas = idxMeas + 2*nMeas;
      end
    end
    
  end %End Lambda Loop

  
  

otherwise
    error(['Unknown forward method: ' pmiModel.Method.Type])
end

%%
%%  Calculate the incident to purturbation ratio
%%
if Debug > 0
   if isfield(pmiModel,'PhiScat')
      pmiModel.IPR = abs(pmiModel.PhiInc) ./ abs(pmiModel.PhiScat);
   end
end
