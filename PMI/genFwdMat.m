%GenFwdMat      Generate the forward matrix from a PMI Model data structure.
%
%   pmiModel = genFwdMat(pmiModel, Debug);
%
%   pmiModel    The Model data structure to updated.
%
%   Debug       The Debug level to execute.
%
%   GenFwdMat generates the forward matrix and incident field from the
%   parameters present in the slab image data structure.  The results are
%   placed in the data structure as the fields A and PhiInc
%   respectively.  Note that this routine automatically moves the diffuse
%   source positions one mean free path into the medium from the positions
%   supplied in pmiModel.
%
%   Calls: none.
%
%   Bugs: assumes diffuse medium exists in the -z domain.
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
%  $Date: 2001/05/24 16:47:09 $
%
%  $Revision: 1.8 $
%
%  $Log: genFwdMat.m,v $
%  Revision 1.8  2001/05/24 16:47:09  dboas
%  The number of coupling coefficients was being counted multiple times when
%  reconstructing musp and mua.  This was fixed.
%
%  Revision 1.7  2001/02/07 16:26:28  dboas
%  Added forward model method type 'ExtBornN'.  This does an iterative extended
%  Born approximation.
%  DON'T FORGET THAT 'BornN' AND 'ExtBornN' DON'T NECESSARILY CONVERGE IF THE
%  PERTURBATION IS TOO LARGE.
%
%  Revision 1.6  2001/02/02 17:51:00  dboas
%  Added Extended Born method.
%
%  Revision 1.5  2000/11/08 18:06:29  dboas
%  Fixed nData to properly handle multiple modulation frequencies.
%
%  Revision 1.4  2000/09/01 19:23:53  dboas
%  Added Model.Method.Type 'FullBorn' and 'BornN'
%  and cleaned up some of the S & D amplitude calculations.
%
%  Revision 1.3  2000/07/27 14:53:41  dboas
%  Added the P(idxLambda) sub-structure for different wavelengths.
%  The optode amplitude must now be a 2 dimensional matrix, where
%     the 1st dimension is optode number and the 2nd is wavelength.
%
%  Revision 1.2  2000/06/27 14:15:34  dboas
%  Small modification to use the correct MeasList
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.5  2000/02/12 00:11:09  dboas
%  The forward matrix now calculates for both the Born-1 and Rytov
%  approximations.
%
%  Revision 3.4  1999/11/18 18:08:07  dboas
%  Fixed an error in how the source and detector were
%  referenced when generating the MeasList
%
%  Revision 3.3  1999/11/18 17:20:39  tgaudett
%  Added Slab
%
%  Revision 3.2  1999/11/16 22:58:55  rjg
%  Added section for slab forward problem, from David.
%
%  Revision 3.1  1999/10/19 21:09:02  rjg
%  Better (although not complete) handling of MeasList.  This is passed on
%  to the DPDWBorn[N/Z]B routines.
%
%  Now requires Lambda data item to be in pmiModel.
%
%  Revision 3.0  1999/06/17 19:29:38  rjg
%  Initial Revision for PMI 3.0
%
%  Revision 2.1  1998/08/20 16:22:01  rjg
%  Moved reseting of SVD flags inside function.
%
%  Revision 2.0  1998/08/05 17:58:36  rjg
%  This function handles calculation of the forward matrice(s) and
%  incident fields.  Moved into a function to make cleanup automatic
%  and handle multiple wavelengths.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmiModel = genFwdMat(pmiModel, Debug)

if nargin < 2
    Debug = 0;
end

%%
%%  Error checking
%%
if ~strcmp(pmiModel.Method.Type, 'Born') & ...
   ~strcmp(pmiModel.Method.Type, 'Rytov') & ...
   ~strcmp(pmiModel.Method.Type, 'FullBorn') & ...
   ~strcmp(pmiModel.Method.Type, 'BornN') & ...
   ~strcmp(pmiModel.Method.Type, 'ExtBorn') & ...
   ~strcmp(pmiModel.Method.Type, 'ExtBornN')
    error(['this algorithm currently only handles a Born, Rytov,' ...
	   ' FullBorn and BornN methods']);
end

%%
%%  Preallocate the memory for the forward matrices and incident field
%%
if strcmp('uniform', lower(pmiModel.CompVol.Type))
    nX = length(pmiModel.CompVol.X);
    nY = length(pmiModel.CompVol.Y);
    nZ = length(pmiModel.CompVol.Z);
    nCompElem = nX * nY *nZ;
else
    nCompElem = size(pmiModel.CompVol.Pos);
end
nVoxs = nCompElem;

nCC = 0;
if isfield(pmiModel.Method,'ObjVec_sd')
  if pmiModel.Method.ObjVec_sd==1
%    nCC = size(pmiModel.Src.Pos,1) + size(pmiModel.Det.Pos,1);
    nCC = size(pmiModel.Src.Pos,1) + size(pmiModel.Det.Pos,1) - 1;
    if pmiModel.ModFreq > 0
      nCC = nCC * 2;
    end
%    nCompElem = nCompElem + nCC;
  end
end

nLambda = length(pmiModel.Lambda);
nFreq = length(pmiModel.ModFreq);
nSrc = size(pmiModel.Src.Pos,1);
nDet = size(pmiModel.Det.Pos,1);

%if ~isfield(pmiModel, 'MeasList')
%    warning('MeasList data item not found, running all permutations');
%    [pSrc nSrc] = getOptodePos(pmiModel.Src);
%    [pDet nDet] = getOptodePos(pmiModel.Det);    
%    pmiModel.MeasList = FullMeasList(nSrc, nDet, nFreq, nLambda);
%end

nData = size(pmiModel.MeasList,1)/(nLambda*nFreq) * ...
    (nFreq + sum(pmiModel.ModFreq > 0));


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
%pmiModel.A = zeros(nData, nCompElem, nLambda);

%%
%%  Loop over wavelength
%%
for idxLambda = 1:nLambda
  if calc_mua & calc_musp
    if strcmpi(pmiModel.Method.Type,'Born') | ...
       strcmpi(pmiModel.Method.Type,'Rytov')
      pmiModel.P(idxLambda).PhiInc = zeros(nData,1);
      pmiModel.P(idxLambda).A = zeros(nData, 2*nCompElem+nCC);
    else
      pmiModel.P(idxLambda).PhiInc_VoxDet = zeros(nCompElem+nCC+nDet,nSrc);
      pmiModel.P(idxLambda).A = zeros(nCompElem+nCC+nDet, 4*nCompElem+nCC+nDet);
    end
  elseif calc_musp
    if strcmpi(pmiModel.Method.Type,'Born') | ...
       strcmpi(pmiModel.Method.Type,'Rytov')
      pmiModel.P(idxLambda).PhiInc = zeros(nData,1);
      pmiModel.P(idxLambda).A = zeros(nData, nCompElem+nCC);
    else
      pmiModel.P(idxLambda).PhiInc_VoxDet = zeros(nCompElem+nCC+nDet,nSrc);
      pmiModel.P(idxLambda).A = zeros(nCompElem+nCC+nDet, 3*nCompElem+nCC+nDet);
    end
  else
    if strcmpi(pmiModel.Method.Type,'Born') | ...
       strcmpi(pmiModel.Method.Type,'Rytov')
      pmiModel.P(idxLambda).PhiInc = zeros(nData,1);
      pmiModel.P(idxLambda).A = zeros(nData, nCompElem+nCC);
    else
      pmiModel.P(idxLambda).PhiInc_VoxDet = zeros(nCompElem+nCC+nDet,nSrc);
      pmiModel.P(idxLambda).A = zeros(nCompElem+nCC+nDet, nCompElem+nCC+nDet);
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
    
    switch lower(pmiModel.Boundary.Geometry)
     case { 'semi-infinite', 'semi', 'extrapolated'}
      if Debug
	fprintf(['Executing extrapolated zero boundary' ...
		 ' computation\n']);
      end
      [A PhiInc] = DPDWBorn1ZB(pmiModel, MeasList, Debug);
      
     case {'infinite', 'inf'}
      if Debug
	fprintf(['Executing infinite medium boundary' ...
		 ' computation\n']);
      end
      [A PhiInc] = DPDWBorn1NB(pmiModel, MeasList, Debug);
      
     case {'slab'}
      [A PhiInc] = DPDWBorn1slab(pmiModel, MeasList, Debug);
      
     case {'slabfd'}
      [A PhiInc] = DPDWBorn1slabFD(pmiModel, MeasList );
      
     otherwise
      error(['Unknown boundary condition: ' pmiModel.Boundary.Geometry]);
    end
    
    %%
    %%  Scale the amplitude of the forward matrix and incident field if 
    %%  necessary.
    %%
    %        if any(pmiModel.Src.Amplitude ~= 1) | any(pmiModel.Det.Amplitude ~= 1)
    %            disp('NEED TO VALIDATE MATRIX AMPLITUDE WEIGHTING!');
    if strcmpi(pmiModel.Method.Type,'Born') | ...
       strcmpi(pmiModel.Method.Type,'Rytov')
      clear SDWeight;
      nMeas = size(MeasList,1);
      for idx = 1:nMeas
	SDWeight(idx,1) = ...
	    pmiModel.Src.Amplitude(MeasList(idx,1),idxLambda) ...
	    * pmiModel.Det.Amplitude(MeasList(idx,2),idxLambda);
      end
      PhiInc = SDWeight .* PhiInc;
    
      if strcmpi(pmiModel.Method.Type,'Born')
	A = rowscale(A, SDWeight);
      end
      
    else

      nSrcs = size(pmiModel.Src.Amplitude,1);
      nDets = size(pmiModel.Det.Amplitude,1);
      nOpt = size(PhiInc,2);
      if nOpt==nSrcs
	PhiInc = rowscale(PhiInc', pmiModel.Src.Amplitude(:, ...
							  idxLambda))';
      else
	% added for Non-linear updates
	PhiInc = rowscale(PhiInc', ...
			  [pmiModel.Src.Amplitude(:, idxLambda);
			   pmiModel.Det.Amplitude(:, idxLambda)]  )';
      end
      
      [nx ny] = size(A);
      nDets = size(pmiModel.Det.Pos,1);
      PhiInc(nx-nDets+1:nx,:) = rowscale(PhiInc(nx-nDets+1:nx,:), ...
				    pmiModel.Det.Amplitude(:,idxLambda) );
      A(nx-nDets+1:nx,:) = rowscale(A(nx-nDets+1:nx,:), ...
				    pmiModel.Det.Amplitude(:,idxLambda) ...
				    );
      
    end
    
      
    
    %       end
    
    %%
    %% If the pmiModel.Method.ObjVec_sd flag is true then add the source
    %% detector amplitude coefficients to the A matrix
    %%
    if isfield(pmiModel.Method,'ObjVec_sd')
      if pmiModel.Method.ObjVec_sd==1
	clear B;
	nMeas = size(MeasList,1);
	nSrcs = size(pmiModel.Src.Pos,1);
	nDets = size(pmiModel.Det.Pos,1);
%	nCC = nSrcs + nDets;
	nCC = nSrcs + nDets - 1; %We assume that the
				 %last detector coupling 
				 %coefficient is known
	if pmiModel.ModFreq > 0 & strcmpi(pmiModel.Method.Type,'Rytov')
	  nCC = nCC * 2;
	end
	B = zeros(nMeas,nCC);
	for idx=1:nMeas
	  if strcmp(pmiModel.Method.Type,'Rytov')
	    B(idx,MeasList(idx,1)) = 1;
	    if pmiModel.ModFreq > 0
%	      B(idx,nSrcs+nDets+MeasList(idx,1)) = sqrt(-1);
	      B(idx,nSrcs+nDets-1+MeasList(idx,1)) = sqrt(-1);
	    end
	  else
	    B(idx,MeasList(idx,1)) = PhiInc(idx)/pmiModel.Src.Amplitude(MeasList(idx,1));
	  end
	  if MeasList(idx,2) ~= nDets
	    if strcmp(pmiModel.Method.Type,'Rytov')
	      B(idx,nSrcs+MeasList(idx,2)) = 1;
	      if pmiModel.ModFreq > 0
%		B(idx,2*nSrcs+nDets+MeasList(idx,2)) = sqrt(-1);
		B(idx,2*nSrcs+nDets-1+MeasList(idx,2)) = sqrt(-1);
	      end
	    else
	      B(idx,nSrcs+MeasList(idx,2)) = PhiInc(idx)/pmiModel.Det.Amplitude(MeasList(idx,2));
	    end
	  end
	end
	
	A = [A B];
      end
    end

    %%
    %%  Copy the forward matrix into the correct block
    %%
    if strcmpi(pmiModel.Method.Type,'Born') | ...
       strcmpi(pmiModel.Method.Type,'Rytov')
      nMeas = size(MeasList, 1);
      pmiModel.P(idxLambda).A(idxRow:idxRow+nMeas-1,:) = real(A);
      pmiModel.P(idxLambda).PhiInc(idxRow:idxRow+nMeas-1,:) = real(PhiInc);
      idxRow = idxRow + nMeas;            
      if pmiModel.ModFreq(idxFreq) > 0
	pmiModel.P(idxLambda).A(idxRow:idxRow+nMeas-1,:) = imag(A);
	pmiModel.P(idxLambda).PhiInc(idxRow:idxRow+nMeas-1,:) = imag(PhiInc);
	idxRow = idxRow + nMeas;
      end
    else
      nMeas = size(A,1);

      % added for non-linear
      pmiModel.P(idxLambda).PhiInc_VoxDet = zeros(nCompElem+nDet,size(PhiInc,2));

      pmiModel.P(idxLambda).A(idxRow:idxRow+nMeas-1,:) = A;
      pmiModel.P(idxLambda).PhiInc_VoxDet(idxRow:idxRow+nMeas-1,:) = PhiInc;
      idxRow = idxRow + nMeas;            
    end

  
  
  end  %END FREQUENCY LOOP

  % MAKE THE MATRIX SPARSE IF 'ExtBornN'
  if strcmpi(pmiModel.Method.Type,'ExtBornN')
    pmiModel.P(idxLambda).A = sparse(pmiModel.P(idxLambda).A);
  end

end   %END WAVELENGTH LOOP

