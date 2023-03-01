%HLM3PTBORN1slab  1st Born approx. of a perturbed Helmholtz equation, z boundary.
%
%   [A Phi_Inc] = Hlm3ptBorn1slab(pmiModel, k, MeasList, zBnd, Debug)
%
%   A           The forward matrix relating the contribution from each voxel to
%               a specific source-detector pair.  Each row is for a different
%               source detector pair.
%
%   Phi_Inc     The incident response at each detector from each source in a
%               column vector.  The same combination pattern as the columns of A.
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
%   HLM3PTBORN1SLAB computes the forward weighting matrix associated with
%   the Born-1 approximation to a spatial varying k.  A zero fluence boudary
%   condition is implemented at the distance zBnd by mirroring the sources
%   and the pertubation responses.  The sources are assumed to be unit
%   amplitude point sources.
%
%   Do not specify the sources or detector exactly at a sampling point to
%   prevent dived by zero errors when evaluating the Green's functions.

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
%  $Date: 2001/08/16 23:03:17 $
%
%  $Revision: 1.11 $
%
%  $Log: Hlm3ptBorn1slab.m,v $
%  Revision 1.11  2001/08/16 23:03:17  dboas
%  Gscat_musp_z had a small mistake in the image.
%
%  Revision 1.10  2001/05/10 22:26:07  dboas
%  The scattered wave from the musp perturbation had the wrong sign
%
%  Revision 1.9  2001/02/02 17:52:37  dboas
%  Removed the 1/3 factor from the calculation of the perturbation matrix for
%  delta_D.  I don't recall
%  why this 1/3 appeared originally, but it shouldn't be there.
%
%  Revision 1.8  2000/11/07 16:47:30  dboas
%  Added v/D scaling for mua perturbation.
%  HAVE NOT YET added SCATTERING perturbation.
%
%  Revision 1.7  2000/09/06 12:07:48  dboas
%  Fixed a small error in the use of multiple images.
%
%  Revision 1.6  2000/08/10 19:03:06  dboas
%  The distances from some of the image voxels to the detector were not being
%  calculated correctly.  This has been fixed.
%
%  Revision 1.5  2000/08/04 16:44:45  tgaudett
%  Fixed extra images
%
%  Revision 1.4  2000/08/04 18:12:07  dboas
%  Using 7 images, updated from 5.
%
%  Revision 1.3  2000/07/27 15:07:59  dboas
%  It is not proper to scale the matrix by (v/D).^2, this is done in
%  DPDWBorn1slab.
%
%  Revision 1.2  2000/06/27 14:24:53  dboas
%  Increased the number of image sources and voxels. Probably not necessary but
%  safe.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.5  2000/01/10 00:14:14  dboas
%  Storing the source and detector lists for use by other functions
%
%  Revision 1.4  1999/11/18 17:50:01  dboas
%  Fixed an error in how the source and detector were
%  referenced in the MeasList
%
%  Revision 1.3  1999/11/18 14:19:01  tgaudett
%  Creation of Slab Geometry Calculation
%
%  Revision 1.2  1999/11/18 13:57:45  tgaudett
%  Creation of Slab Geometry Calculation.
%  With the correct log file info
%
% 
%  Revision 1.1  1998/06/03 16:04:48  DAB
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = Hlm3ptBorn1slab(pmiModel, k, MeasList, zBnd, Debug)

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
  A = zeros(2*nPts, nMeas) + j * ones(2*nPts, nMeas);
else
  A = zeros(nPts, nMeas) + j * ones(nPts, nMeas);
end
Phi_Inc = zeros(nMeas, 1) + j * ones(nMeas, 1);

D = pmiModel.v(idxLambda) / (3 * pmiModel.Mu_sp(idxLambda));

%%
%%  Loop over each source detector combination
%%
if Debug
    fprintf('Meas #: ');
end


for iMeas = 1:nMeas
    if Debug
        fprintf('%d  ', iMeas);
    end
    Src = pSrc(MeasList(iMeas, 1),:);
    Det = pDet(MeasList(iMeas, 2),:);

    %%
    %%  Compute the incident response on every voxel from each measurement 
    %%  (source detector pair).  The source is represented by a point source.
    %%
    rho_sq = (Xm - Src(1)).^2 + (Ym - Src(2)).^2;
    rs = sqrt(rho_sq + (Zm - Src(3)).^2);
    
    %%
    %%  Add in the image source due to the z=0 boundary
    %%
    z1i = (-1) * sign(Thickness) * (2 * zBnd + abs(Src(3)));
    rs1i = sqrt(rho_sq + (Zm - z1i).^2);
    
    %%
    %%  Add in the image source due to the z=THICKNESS boundary
    %%
    z2i = (2*Thickness + sign(Thickness)* 2 * zBnd - Src(3));
    rs2i = sqrt(rho_sq + (Zm - z2i).^2);
    
    %%
    %% Image of z=0 Image source due to z=THICKNESS boundary
    %%
    z1ii = (2*Thickness + sign(Thickness)* 2 * zBnd - z1i);
    rs1ii = sqrt(rho_sq + (Zm - z1ii).^2);
    
    % Next Images
    z2ii = (-sign(Thickness) * 2 * zBnd - z2i);
    rs2ii = sqrt(rho_sq + (Zm - z2ii).^2);
      
    z1iii = (-sign(Thickness) * 2 * zBnd - z1ii);
    rs1iii = sqrt(rho_sq + (Zm - z1iii).^2);
    
    %%
    %%  Next Images
    %%
    z2iii = (2*Thickness + sign(Thickness)* 2 * zBnd - z2ii);
    rs2iii = sqrt(rho_sq + (Zm - z2iii).^2);
    
    z1iiii = (2*Thickness + sign(Thickness)* 2 * zBnd - z1iii);
    rs1iiii = sqrt(rho_sq + (Zm - z1iiii).^2);
    
    if calc_mua
      phi_inc_mua = exp(j * k * rs) ./ (-4 * pi * rs) ...
	  - exp(j * k * rs1i) ./ (-4 * pi * rs1i) ...
	  - exp(j * k * rs2i) ./ (-4 * pi * rs2i) ...
	  + exp(j * k * rs1ii) ./ (-4 * pi * rs1ii) ...
	  + exp(j * k * rs2ii) ./ (-4 * pi * rs2ii) ...
	  - exp(j * k * rs1iii) ./ (-4 * pi * rs1iii)...
	  - exp(j * k * rs2iii) ./ (-4 * pi * rs2iii) ...
	  + exp(j * k * rs1iiii) ./ (-4 * pi * rs1iiii);
    end
    if calc_musp
      phi_inc_musp_x = (j*k -1./rs).*(Xm-Src(1)).*exp(j * k * rs) ./ (-4 * pi * rs.^2) ...
	  - (j*k -1./rs1i).*(Xm-Src(1)).*exp(j * k * rs1i) ./ (-4 * pi * rs1i.^2) ...
	  - (j*k -1./rs2i).*(Xm-Src(1)).*exp(j * k * rs2i) ./ (-4 * pi * rs2i.^2) ...
	  + (j*k -1./rs1ii).*(Xm-Src(1)).*exp(j * k * rs1ii) ./ (-4 * pi * rs1ii.^2)...
	  + (j*k -1./rs2ii).*(Xm-Src(1)).*exp(j * k * rs2ii) ./ (-4 * pi * rs2ii.^2) ...
	  - (j*k -1./rs1iii).*(Xm-Src(1)).*exp(j * k * rs1iii) ./ (-4 * pi * rs1iii.^2)...
	  - (j*k -1./rs2iii).*(Xm-Src(1)).*exp(j * k * rs2iii) ./ (-4 * pi * rs2iii.^2) ...
	  + (j*k -1./rs1iiii).*(Xm-Src(1)).*exp(j * k * rs1iiii) ./ (-4 * pi * rs1iiii.^2);
    
      phi_inc_musp_y = (j*k -1./rs).*(Ym-Src(2)).*exp(j * k * rs) ./ (-4 * pi * rs.^2) ...
	  - (j*k -1./rs1i).*(Ym-Src(2)).*exp(j * k * rs1i) ./ (-4 * pi * rs1i.^2) ...
	  - (j*k -1./rs2i).*(Ym-Src(2)).*exp(j * k * rs2i) ./ (-4 * pi * rs2i.^2) ...
	  + (j*k -1./rs1ii).*(Ym-Src(2)).*exp(j * k * rs1ii) ./ (-4 * pi * rs1ii.^2)...
	  + (j*k -1./rs2ii).*(Ym-Src(2)).*exp(j * k * rs2ii) ./ (-4 * pi * rs2ii.^2) ...
	  - (j*k -1./rs1iii).*(Ym-Src(2)).*exp(j * k * rs1iii) ./ (-4 * pi * rs1iii.^2)...
	  - (j*k -1./rs2iii).*(Ym-Src(2)).*exp(j * k * rs2iii) ./ (-4 * pi * rs2iii.^2) ...
	  + (j*k -1./rs1iiii).*(Ym-Src(2)).*exp(j * k * rs1iiii) ./ (-4 * pi * rs1iiii.^2);
    
      phi_inc_musp_z = (j*k -1./rs).*(Zm-Src(3)).*exp(j * k * rs) ./ (-4 * pi * rs.^2) ...
	  - (j*k -1./rs1i).*(Zm-z1i).*exp(j * k * rs1i) ./ (-4 * pi * rs1i.^2) ...
	  - (j*k -1./rs2i).*(Zm-z2i).*exp(j * k * rs2i) ./ (-4 * pi * rs2i.^2) ...
	  + (j*k -1./rs1ii).*(Zm-z1ii).*exp(j * k * rs1ii) ./ (-4 * pi * rs1ii.^2)...
	  + (j*k -1./rs2ii).*(Zm-z2ii).*exp(j * k * rs2ii) ./ (-4 * pi * rs2ii.^2) ...
	  - (j*k -1./rs1iii).*(Zm-z1iii).*exp(j * k * rs1iii) ./ (-4 * pi * rs1iii.^2)...
	  - (j*k -1./rs2iii).*(Zm-z2iii).*exp(j * k * rs2iii) ./ (-4 * pi * rs2iii.^2) ...
	  + (j*k -1./rs1iiii).*(Zm-z1iiii).*exp(j * k * rs1iiii) ./ (-4 * pi * rs1iiii.^2);
    end
    
    
    %%
    %%  Compute the incident response at this detector  from all
    %%  of the voxels.  Note that the negative sign in the Green's function
    %%  is cancelled by the negative sign in the Born approximation
    %%
    rho_sq = (Xm - Det(1)).^2 + (Ym - Det(2)).^2;
    rd = sqrt(rho_sq + (Zm - Det(3)).^2);

    %%
    %%  Add in the image of the voxel due to the z=0 boundary
    %%
    z1i = (-1)*sign(Thickness) * (2 * zBnd + abs(Zm));
    rd1i = sqrt(rho_sq + (Det(3) - z1i).^2);
    
    %%
    %%  Add in the image of the voxel due to the z=THICKNESS boundary
    %%
    z2i = (2*Thickness + sign(Thickness)* 2 * zBnd - Zm);
    rd2i = sqrt(rho_sq + (Det(3) - z2i).^2);
    
    %%
    %% Image of z=0 Image voxel due to z=THICKNESS boundary
    %%
    z1ii = (2*Thickness + sign(Thickness)* 2 * zBnd - z1i);
    rd1ii = sqrt(rho_sq + (Det(3) - z1ii).^2);
    
    %%
    %% Next images
    %%
    z2ii = (-sign(Thickness) * 2 * zBnd - z2i);
    rd2ii = sqrt(rho_sq + (Det(3) - z2ii).^2 );
    
    z1iii = (-sign(Thickness) * 2 * zBnd - z1ii);
    rd1iii = sqrt(rho_sq + (Det(3) - z1iii).^2 );
    
    %%
    %%  Next Images
    %%
    z2iii = (2*Thickness + sign(Thickness)* 2 * zBnd - z2ii);
    rd2iii = sqrt(rho_sq + (Det(3) - z2iii).^2);
    
    z1iiii = (2*Thickness + sign(Thickness)* 2 * zBnd - z1iii);
    rd1iiii = sqrt(rho_sq + (Det(3) - z1iiii).^2);
    

    if calc_mua
      Gscat_mua = (exp(j * k * rd) ./ (4 * pi * rd) ...
		   - exp(j * k * rd1i) ./ (4 * pi * rd1i) ...
		   - exp(j * k * rd2i) ./ (4 * pi * rd2i) ...
		   + exp(j * k * rd1ii) ./ (4 * pi * rd1ii) ...
		   + exp(j * k * rd2ii) ./ (4 * pi * rd2ii) ...
		   - exp(j * k * rd1iii) ./ (4 * pi * rd1iii) ...
		   - exp(j * k * rd2iii) ./ (4 * pi * rd2iii) ...
		   + exp(j * k * rd1iiii) ./ (4 * pi * rd1iiii) ...
		   ) * volVoxel * pmiModel.v(idxLambda) / D;
    
      A(1:nPts, iMeas) = Gscat_mua .* phi_inc_mua;
    end    

    if calc_musp
      Gscat_musp_x = ((j*k -1./rd).*(Xm-Det(1)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
	  - (j*k -1./rd1i).*(Xm-Det(1)).*exp(j * k * rd1i) ./ (4 * pi * rd1i.^2) ...
	  - (j*k -1./rd2i).*(Xm-Det(1)).*exp(j * k * rd2i) ./ (4 * pi * rd2i.^2) ...
	  + (j*k -1./rd1ii).*(Xm-Det(1)).*exp(j * k * rd1ii) ./ (4 * pi * rd1ii.^2)...
	  + (j*k -1./rd2ii).*(Xm-Det(1)).*exp(j * k * rd2ii) ./ (4 * pi * rd2ii.^2) ...
	  - (j*k -1./rd1iii).*(Xm-Det(1)).*exp(j * k * rd1iii) ./ (4 * pi * rd1iii.^2)...
	  - (j*k -1./rd2iii).*(Xm-Det(1)).*exp(j * k * rd2iii) ./ (4 * pi * rd2iii.^2) ...
	  + (j*k -1./rd1iiii).*(Xm-Det(1)).*exp(j * k * rd1iiii) ./ (4 * pi * rd1iiii.^2)...
		      ) * volVoxel;
      
      Gscat_musp_y = ((j*k -1./rd).*(Ym-Det(2)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
	  - (j*k -1./rd1i).*(Ym-Det(2)).*exp(j * k * rd1i) ./ (4 * pi * rd1i.^2) ...
	  - (j*k -1./rd2i).*(Ym-Det(2)).*exp(j * k * rd2i) ./ (4 * pi * rd2i.^2) ...
	  + (j*k -1./rd1ii).*(Ym-Det(2)).*exp(j * k * rd1ii) ./ (4 * pi * rd1ii.^2)...
	  + (j*k -1./rd2ii).*(Ym-Det(2)).*exp(j * k * rd2ii) ./ (4 * pi * rd2ii.^2) ...
	  - (j*k -1./rd1iii).*(Ym-Det(2)).*exp(j * k * rd1iii) ./ (4 * pi * rd1iii.^2)...
	  - (j*k -1./rd2iii).*(Ym-Det(2)).*exp(j * k * rd2iii) ./ (4 * pi * rd2iii.^2) ...
	  + (j*k -1./rd1iiii).*(Ym-Det(2)).*exp(j * k * rd1iiii) ./ (4 * pi * rd1iiii.^2)...
		      ) * volVoxel;

      Gscat_musp_z = ((j*k -1./rd).*(Zm-Det(3)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
		      - (j*k -1./rd1i).*(z1i-Det(3)).*exp(j * k * rd1i) ./ (4 * pi * rd1i.^2) ...
		      - (j*k -1./rd2i).*(z2i-Det(3)).*exp(j * k * rd2i) ./ (4 * pi * rd2i.^2) ...
		      + (j*k -1./rd1ii).*(z1ii-Det(3)).*exp(j * k * rd1ii) ./ (4 * pi * rd1ii.^2)...
		      + (j*k -1./rd2ii).*(z2ii-Det(3)).*exp(j * k * rd2ii) ./ (4 * pi * rd2ii.^2) ...
		      - (j*k -1./rd1iii).*(z1iii-Det(3)).*exp(j * k * rd1iii) ./ (4 * pi * rd1iii.^2)...
		      - (j*k -1./rd2iii).*(z2iii-Det(3)).*exp(j * k * rd2iii) ./ (4 * pi * rd2iii.^2) ...
		      + (j*k -1./rd1iiii).*(z1iiii-Det(3)).*exp(j * k * rd1iiii) ./ (4 * pi * rd1iiii.^2)...
		      ) * volVoxel;

    
      if ~calc_mua
	A(1:nPts, iMeas) = (Gscat_musp_x .* phi_inc_musp_x ...
	    + Gscat_musp_y .* phi_inc_musp_y ...
	    + Gscat_musp_z .* phi_inc_musp_z);
      else
	A(nPts+1:2*nPts, iMeas) = (Gscat_musp_x .* phi_inc_musp_x ...
	    + Gscat_musp_y .* phi_inc_musp_y ...
	    + Gscat_musp_z .* phi_inc_musp_z);
      end

    end

    %%
    %%  Compute the incident response at the detector.
    %%
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
    rsrcdet = norm(Src - Det);
    rimg1det = norm(Image1Src - Det);
    rimg2det = norm(Image2Src - Det);
    rimg1iidet = norm(Image1iiSrc - Det);
    rimg2iidet = norm(Image2iiSrc - Det);
    rimg1iiidet = norm(Image1iiiSrc - Det);
    rimg2iiidet = norm(Image2iiiSrc - Det);
    rimg1iiiidet = norm(Image1iiiiSrc - Det);
    Phi_Inc(iMeas) = exp(j * k * rsrcdet) ./ (-4 * pi * rsrcdet) - ...
       exp(j * k * rimg1det) ./ (-4 * pi * rimg1det) - ...
       exp(j * k * rimg2det) ./ (-4 * pi * rimg2det) + ...
       exp(j * k * rimg1iidet) ./ (-4 * pi * rimg1iidet) + ...
       exp(j * k * rimg2iidet) ./ (-4 * pi * rimg2iidet) - ...
       exp(j * k * rimg1iiidet) ./ (-4 * pi * rimg1iiidet) - ...
       exp(j * k * rimg2iiidet) ./ (-4 * pi * rimg2iiidet) + ...
       exp(j * k * rimg1iiiidet) ./ (-4 * pi * rimg1iiiidet);

end

A = A.';
if Debug
    fprintf('\n');
end
