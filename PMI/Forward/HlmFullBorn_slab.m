%HLMFullBorn_slab
%
%   Calculates the perturbation matrix for the Nth Born Approximation
%
%   [A Phi_Inc] = Hlm3ptBorn1slab(pmiModel, k, MeasList, zBnd, Debug)
%
%   A           The forward matrix relating the contribution from each voxel to
%               a specific voxel or detector.  The matrix is
%               square.  Each col/row is for a different
%               voxel or detector pair.
%
%   Phi_Inc     The incident response at each voxel or detector
%               from each source.  Each source is in a different
%               column.  That is, the matrix is #voxels + #dets by
%               #srcs.  Each row is
%               the same as the cols/rows of A.
%
%   pmiModel    The PMI Model structure contain the following fields: CompVol,
%               Mu_sp, Mu_a, v, idxRefr, and f.
%
%   k           The complex wavenumber for the Helmholtz equation.
%
%   MeasList    The measurement list for this simulation.  All frequencies and
%               wavelengths must be the same.
%
%   zBnd        The position of the zero value boundary (Dirichlet
%               conditions).
%
%   Debug       OPTIONAL: Print out debugging info.
%
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
%  $Date: 2001/02/20 16:31:50 $
%
%  $Revision: 1.5 $
%
%  $Log: HlmFullBorn_slab.m,v $
%  Revision 1.5  2001/02/20 16:31:50  dboas
%  Fixed a typo for checking 'ExtBorn' method.
%
%  Revision 1.4  2001/02/07 16:30:41  dboas
%  Added forward model method type 'ExtBornN'.  This does an iterative extended
%  Born approximation.
%  DON'T FORGET THAT 'BornN' AND 'ExtBornN' DON'T NECESSARILY CONVERGE IF THE
%  PERTURBATION IS TOO LARGE.
%
%  Removed 1/3 factor from the matrix.
%
%  Revision 1.3  2001/02/02 17:52:57  dboas
%  Added extended born with radius of influence which is specified by
%     pmiModel.Method.ExtBorn_Radius
%  NOTE THAT EXTENDED BORN AND FULL BORN PRESENTLY DO NOT WORK FOR DELTA_D.
%
%  Revision 1.2  2000/11/07 16:50:22  dboas
%  Scale mua perturbation by v/D.
%  Deal with self-voxel contribution for absorption perturbation.
%  Not yet dealing with self-voxel for scattering perturbation since I am not
%  convinced that it exists.
%
%  Revision 1.1  2000/09/06 12:17:40  dboas
%  Initial version.  Calculates the perturbation matrix for the Nth
%  Born Approximation.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = HlmFullBorn_slab(pmiModel, k, MeasList, zBnd, Debug)

if nargin < 5
    Debug = 0;
end
nMeas = size(MeasList,1);
idxLambda = MeasList(1,4);

D = pmiModel.v(idxLambda) / (3 * pmiModel.Mu_sp(idxLambda));

%%
%% Extract the slab thickness
%%
Thickness = pmiModel.Boundary.Thickness;

%%
%%  Extract the source and detector positions
%%  Move the effective source position 1 mean free path into the medium
%%
pSrc = getOptodePos(pmiModel.Src);
foo = find(pSrc(:,3)==0);
pSrc(foo,3) = pSrc(foo,3) + sign(Thickness) * (1/pmiModel.Mu_sp(idxLambda));
foo = find(pSrc(:,3)==Thickness);
pSrc(foo,3) = pSrc(foo,3) - sign(Thickness) * (1/pmiModel.Mu_sp(idxLambda));

pDet = getOptodePos(pmiModel.Det);
nDets = size(pmiModel.Det.Pos,1);
nSrcs = size(pmiModel.Src.Pos,1);

%%
%%  Create the sampling volume
%%
if strcmp(pmiModel.CompVol.Type, 'uniform') == 0
    error('This routine only handles uniform voxelations');
end
[Xm Ym Zm] = meshgrid(pmiModel.CompVol.X, pmiModel.CompVol.Y, pmiModel.CompVol.Z);
Xm = Xm(:);
Ym = Ym(:);
Zm = Zm(:);
nPts = length(Xm);
volVoxel =  pmiModel.CompVol.XStep * pmiModel.CompVol.YStep * ...
    pmiModel.CompVol.ZStep;
R = (volVoxel / (4*pi/3))^(1/3);

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
%%  A is created transposed for better caching performance
%%
if calc_mua & calc_musp
  num1 = 4*nPts+nDets;
  num2 = nPts+nDets;
  A = zeros(num1, num2) + j * zeros(num1, num2);
elseif calc_musp
  num1 = 3*nPts+nDets;
  num2 = nPts+nDets;
  A = zeros(num1, num2) + j * zeros(num1, num2);
else
  num1 = nPts+nDets;
  num2 = nPts+nDets;
  A = zeros(num1, num2) + j * zeros(num1, num2);
end
Phi_Inc = zeros(num2, nSrcs) + j * ones(num2, nSrcs);

%%
%%  Loop over each voxel and detector
%%
if Debug
    fprintf('Meas #: ');
end


for idx = 1:num2
  if Debug
    fprintf('%d  ', idx);
  end
  
  if idx <= (num2-nDets)
    Pos = [Xm(idx) Ym(idx) Zm(idx)];
  else
    Pos = pDet(idx-(num2-nDets),:);
  end
  
  
  %%
  %%  Compute the Green's function from each voxel to each voxel/detector.
  %%  Note that the negative sign in the Green's function
  %%  is cancelled by the negative sign in the Born approximation
  %%
  rho_sq = (Xm - Pos(1)).^2 + (Ym - Pos(2)).^2;
  rd = sqrt(rho_sq + (Zm - Pos(3)).^2);
  
  %%
  %%  Add in the image of the voxel due to the z=0 boundary
  %%
  z1i = (-1)*sign(Thickness) * (2 * zBnd + abs(Zm));
  rd1i = sqrt(rho_sq + (Pos(3) - z1i).^2);
  
  %%
  %%  Add in the image of the voxel due to the z=THICKNESS boundary
  %%
  z2i = (2*Thickness + sign(Thickness)* 2 * zBnd - Zm);
  rd2i = sqrt(rho_sq + (Pos(3) - z2i).^2);
  
  %%
  %% Image of z=0 Image voxel due to z=THICKNESS boundary
  %%
  z1ii = (2*Thickness + sign(Thickness)* 2 * zBnd - z1i);
  rd1ii = sqrt(rho_sq + (Pos(3) - z1ii).^2);
  
  %%
  %% Next images
  %%
  z2ii = (-sign(Thickness) * 2 * zBnd - z2i);
  rd2ii = sqrt(rho_sq + (Pos(3) - z2ii).^2 );
  
  z1iii = (-sign(Thickness) * 2 * zBnd - z1ii);
  rd1iii = sqrt(rho_sq + (Pos(3) - z1iii).^2 );
  
  %%
  %%  Next Images
  %%
  z2iii = (2*Thickness + sign(Thickness)* 2 * zBnd - z2ii);
  rd2iii = sqrt(rho_sq + (Pos(3) - z2iii).^2);
  
  z1iiii = (2*Thickness + sign(Thickness)* 2 * zBnd - z1iii);
  rd1iiii = sqrt(rho_sq + (Pos(3) - z1iiii).^2);
  

  if ~(strcmpi(pmiModel.Method.Type,'ExtBorn') | ...
       strcmpi(pmiModel.Method.Type,'ExtBornN') ) | idx>(num2-nDets)
    nonZeroIdx = find(rd~=0);
  else
    nonZeroIdx = find(rd~=0 & rd<=pmiModel.Method.ExtBorn_Radius);
  end
  ZeroIdx = find(rd==0);

  if calc_mua
    Gscat_mua = zeros(nPts,1);
    Gscat_mua(nonZeroIdx) = (exp(j * k * rd(nonZeroIdx)) ./ (4 * pi ...
				* rd(nonZeroIdx)) )*volVoxel*... 
	pmiModel.v(idxLambda) / D; 
    Gscat_mua(ZeroIdx) = (3 * pmiModel.Mu_sp(idxLambda)) * ...
	(-1/(k^2)) * (1-(1-j*k*R)*exp(j*k*R))*...
	pmiModel.v(idxLambda) / D; 

    Gscat_mua = Gscat_mua + (...
    		 - exp(j * k * rd1i) ./ (4 * pi * rd1i) ...
		 - exp(j * k * rd2i) ./ (4 * pi * rd2i) ...
		 + exp(j * k * rd1ii) ./ (4 * pi * rd1ii) ...
		 + exp(j * k * rd2ii) ./ (4 * pi * rd2ii) ...
		 - exp(j * k * rd1iii) ./ (4 * pi * rd1iii) ...
		 - exp(j * k * rd2iii) ./ (4 * pi * rd2iii) ...
		 + exp(j * k * rd1iiii) ./ (4 * pi * rd1iiii) ...
		 ) * volVoxel;
    
    A(1:nPts, idx) = Gscat_mua;
  end    
  
  if calc_musp
    Gscat_musp_x = zeros(nPts,1);
    Gscat_musp_y = zeros(nPts,1);
    Gscat_musp_z = zeros(nPts,1);

    Gscat_musp_x(nonZeroIdx) = ((j*k -1./rd(nonZeroIdx)).* ...
	(Xm(nonZeroIdx)-Pos(1)).*exp(j * k * rd(nonZeroIdx)) ./ (-4 ...
			  * pi * rd(nonZeroIdx).^2) )*volVoxel;  
    
    Gscat_musp_x = Gscat_musp_x + (...
		    - (j*k -1./rd1i).*(Xm-Pos(1)).*exp(j * k * rd1i) ./ (-4 * pi * rd1i.^2) ...
		    - (j*k -1./rd2i).*(Xm-Pos(1)).*exp(j * k * rd2i) ./ (-4 * pi * rd2i.^2) ...
		    + (j*k -1./rd1ii).*(Xm-Pos(1)).*exp(j * k * rd1ii) ./ (-4 * pi * rd1ii.^2)...
		    + (j*k -1./rd2ii).*(Xm-Pos(1)).*exp(j * k * rd2ii) ./ (-4 * pi * rd2ii.^2) ...
		    - (j*k -1./rd1iii).*(Xm-Pos(1)).*exp(j * k * rd1iii) ./ (-4 * pi * rd1iii.^2)...
		    - (j*k -1./rd2iii).*(Xm-Pos(1)).*exp(j * k * rd2iii) ./ (-4 * pi * rd2ii.^2) ...
		    + (j*k -1./rd1iiii).*(Xm-Pos(1)).*exp(j * k * rd1iiii) ./ (-4 * pi * rd1iiii.^2)...
		    ) * volVoxel;
    
    Gscat_musp_y(nonZeroIdx) = ((j*k -1./rd(nonZeroIdx)).* ...
	(Ym(nonZeroIdx)-Pos(2)).*exp(j * k * rd(nonZeroIdx)) ./ (-4 ...
			  * pi * rd(nonZeroIdx).^2) )*volVoxel;  
    
    Gscat_musp_y = Gscat_musp_y + (...
		    - (j*k -1./rd1i).*(Ym-Pos(2)).*exp(j * k * rd1i) ./ (-4 * pi * rd1i.^2) ...
		    - (j*k -1./rd2i).*(Ym-Pos(2)).*exp(j * k * rd2i) ./ (-4 * pi * rd2i.^2) ...
		    + (j*k -1./rd1ii).*(Ym-Pos(2)).*exp(j * k * rd1ii) ./ (-4 * pi * rd1ii.^2)...
		    + (j*k -1./rd2ii).*(Ym-Pos(2)).*exp(j * k * rd2ii) ./ (-4 * pi * rd2ii.^2) ...
		    - (j*k -1./rd1iii).*(Ym-Pos(2)).*exp(j * k * rd1iii) ./ (-4 * pi * rd1iii.^2)...
		    - (j*k -1./rd2iii).*(Ym-Pos(2)).*exp(j * k * rd2iii) ./ (-4 * pi * rd2ii.^2) ...
		    + (j*k -1./rd1iiii).*(Ym-Pos(2)).*exp(j * k * rd1iiii) ./ (-4 * pi * rd1iiii.^2)...
		    ) * volVoxel;
    
    Gscat_musp_z(nonZeroIdx) = ((j*k -1./rd(nonZeroIdx)).* ...
	(Zm(nonZeroIdx)-Pos(3)).*exp(j * k * rd(nonZeroIdx)) ./ (-4 ...
			  * pi * rd(nonZeroIdx).^2) )*volVoxel;  
    
    Gscat_musp_z = Gscat_musp_z + (...
		    - (j*k -1./rd1i).*(Zm-z1i).*exp(j * k * rd1i) ./ (-4 * pi * rd1i.^2) ...
		    - (j*k -1./rd2i).*(Zm-z2i).*exp(j * k * rd2i) ./ (-4 * pi * rd2i.^2) ...
		    + (j*k -1./rd1ii).*(Zm-z1ii).*exp(j * k * rd1ii) ./ (-4 * pi * rd1ii.^2)...
		    + (j*k -1./rd2ii).*(Zm-z2ii).*exp(j * k * rd2ii) ./ (-4 * pi * rd2ii.^2) ...
		    - (j*k -1./rd1iii).*(Zm-z1iii).*exp(j * k * rd1iii) ./ (-4 * pi * rd1iii.^2)...
		    - (j*k -1./rd2iii).*(Zm-z2iii).*exp(j * k * rd2iii) ./ (-4 * pi * rd2ii.^2) ...
		    + (j*k -1./rd1iiii).*(Zm-z1iiii).*exp(j * k * rd1iiii) ./ (-4 * pi * rd1iiii.^2)...
		    ) * volVoxel;    
    
    if ~calc_mua
      A(1:nPts, idx) = Gscat_musp_x;
      A(nPts+1:2*nPts, idx) = Gscat_musp_y;
      A(2*nPts+1:3*nPts, idx) = Gscat_musp_z;
    else
      A(nPts+1:2*nPts, idx) = Gscat_musp_x;
      A(2*nPts+1:3*nPts, idx) = Gscat_musp_y;
      A(3*nPts+1:4*nPts, idx) = Gscat_musp_z;
    end

  end


  %%
  %%  Compute the incident response at the voxel/detector from each source.
  %%
  for idxSrc = 1:nSrcs
    Src = pSrc(idxSrc,:);
  
    z1i = -1*sign(Thickness)*(2*zBnd+abs(Src(3)));
    z2i = (2*Thickness+sign(Thickness)*2*zBnd-Src(3));
    z1ii = (2*Thickness + sign(Thickness) * 2 * zBnd - z1i);
    z2ii = -1*sign(Thickness)*(2*zBnd+abs(z2i));
    z1iii = -1*sign(Thickness)*(2*zBnd+abs(z1ii));
    z2iii = (2*Thickness+sign(Thickness)*2*zBnd-z2ii);
    z1iiii = (2*Thickness + sign(Thickness) * 2 * zBnd - z1iii);
    Image1Src = [Src(1) Src(2) z1i];
    Image2Src = [Src(1) Src(2) z2i];
    Image1iiSrc = [Src(1) Src(2) z1ii];
    Image2iiSrc = [Src(1) Src(2) z2ii];
    Image1iiiSrc = [Src(1) Src(2) z1iii];
    Image2iiiSrc = [Src(1) Src(2) z2iii];
    Image1iiiiSrc = [Src(1) Src(2) z1iiii];
    rsrcdet = norm(Src - Pos);
    rimg1det = norm(Image1Src - Pos);
    rimg2det = norm(Image2Src - Pos);
    rimg1iidet = norm(Image1iiSrc - Pos);
    rimg2iidet = norm(Image2iiSrc - Pos);
    rimg1iiidet = norm(Image1iiiSrc - Pos);
    rimg2iiidet = norm(Image2iiiSrc - Pos);
    rimg1iiiidet = norm(Image1iiiiSrc - Pos);
    Phi_Inc(idx,idxSrc) = exp(j * k * rsrcdet) ./ (-4 * pi * rsrcdet) - ...
	exp(j * k * rimg1det) ./ (-4 * pi * rimg1det) - ...
	exp(j * k * rimg2det) ./ (-4 * pi * rimg2det) + ...
	exp(j * k * rimg1iidet) ./ (-4 * pi * rimg1iidet) + ...
	exp(j * k * rimg2iidet) ./ (-4 * pi * rimg2iidet) - ...
	exp(j * k * rimg1iiidet) ./ (-4 * pi * rimg1iiidet) - ...
	exp(j * k * rimg2iiidet) ./ (-4 * pi * rimg2iiidet) + ...
	exp(j * k * rimg1iiiidet) ./ (-4 * pi * rimg1iiiidet);
    
  end
  
end



A = A.';
if Debug
    fprintf('\n');
end
