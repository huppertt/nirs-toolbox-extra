%HLM3PTBORN1ZB  1st Born approx. of a perturbed Helmholtz equation, z boundary.
%
%   [A Phi_Inc] = Hlm3ptBorn1ZB(pmiModel, k, MeasList, zBnd, Debug)
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
%   HLM3PTBORN1ZB computes the forward weighting matrix associated with
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
%  $Revision: 1.6 $
%
%  $Log: Hlm3ptBorn1ZB.m,v $
%  Revision 1.6  2001/08/16 23:03:17  dboas
%  Gscat_musp_z had a small mistake in the image.
%
%  Revision 1.5  2001/05/10 22:26:07  dboas
%  The scattered wave from the musp perturbation had the wrong sign
%
%  Revision 1.4  2001/02/02 17:52:37  dboas
%  Removed the 1/3 factor from the calculation of the perturbation matrix for
%  delta_D.  I don't recall
%  why this 1/3 appeared originally, but it shouldn't be there.
%
%  Revision 1.3  2000/11/07 16:47:30  dboas
%  Added v/D scaling for mua perturbation.
%  HAVE NOT YET added SCATTERING perturbation.
%
%  Revision 1.2  2000/07/27 15:06:42  dboas
%  Added calculation for scattering perturbation.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.2  1999/11/10 20:23:07  rjg
%  Fixed bug in MeasList that swapped source and detector columns (from David).
%
%  Revision 3.1  1999/10/19 21:58:12  rjg
%  Major change in the handling of the measurement sequence.  The sequence is now
%  defined by MeasList and is no longer all detectors and all sources.
%
%  Revision 3.0  1999/06/17 17:39:56  rjg
%  Initial Revision for PMI 3.0
%
%  Revision 2.1  1999/02/05 20:55:57  rjg
%  Clearified help section.
%
%  Revision 2.0  1998/08/05 15:50:55  rjg
%  Start of version 2
%  Handles structure for computational volume.
%
%  Revision 1.3  1998/06/18 14:56:44  rjg
%  Transformed x,y,z locations into column vectors so that it does not
%  have to be done later in functions of these variables.
%  Changed the variable Phi_Det to Phi_Inc since it is the incident field.
%
%  Revision 1.2  1998/06/05 18:12:20  rjg
%  Changed Psi to Phi.
%
%  Revision 1.1  1998/06/03 16:06:42  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = Hlm3ptBorn1ZB(pmiModel, k, MeasList, zBnd, Debug)

if nargin < 5
    Debug = 0;
end
nMeas = size(MeasList,1);
idxLambda = MeasList(1,4);

%%
%%  Extract the source and detector positions
%%  Move the effective source position 1 mean free path into the medium
%%
pSrc = getOptodePos(pmiModel.Src);
pSrc(:,3) = pSrc(:,3)-(1/pmiModel.Mu_sp(idxLambda));
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
    %%  Add in the image source due to the boundary
    %%
    z1i = 2 * zBnd + abs(Src(3));
    rs1i = sqrt(rho_sq + (Zm - z1i).^2);


    if calc_mua
      phi_inc_mua = exp(j * k * rs) ./ (-4 * pi * rs) ...
	  - exp(j * k * rs1i) ./ (-4 * pi * rs1i);
    end
    if calc_musp
      phi_inc_musp_x = (j*k -1./rs).*(Xm-Src(1)).*exp(j * k * rs) ./ (-4 * pi * rs.^2) ...
	  - (j*k -1./rs1i).*(Xm-Src(1)).*exp(j * k * rs1i) ./ (-4 * pi * rs1i.^2);
    
      phi_inc_musp_y = (j*k -1./rs).*(Ym-Src(2)).*exp(j * k * rs) ./ (-4 * pi * rs.^2) ...
	  - (j*k -1./rs1i).*(Ym-Src(2)).*exp(j * k * rs1i) ./ (-4 * pi * rs1i.^2);
    
      phi_inc_musp_z = (j*k -1./rs).*(Zm-Src(3)).*exp(j * k * rs) ./ (-4 * pi * rs.^2) ...
	  - (j*k -1./rs1i).*(Zm-z1i).*exp(j * k * rs1i) ./ (-4 * pi * rs1i.^2);
    end

    %%
    %%  Compute the incident response at this detector  from all
    %%  of the voxels.  Note that the negative sign in the Green's function
    %%  is cancelled by the negative sign in the Born approximation
    %%
    rho_sq = (Xm - Det(1)).^2 + (Ym - Det(2)).^2;
    rd = sqrt(rho_sq + (Zm - Det(3)).^2);
%    Gscat = exp(i * k * r) ./ (4 * pi * r) * volVoxel;

    %%
    %%  Add in the image of the voxel
    %%
    z1i = 2 * zBnd + abs(Zm);
    rd1i = sqrt(rho_sq + (Det(3) - z1i).^2);
%    Gscat = Gscat - exp(i * k * r) ./ (4 * pi * r) * volVoxel;
%    A(:, iMeas) = Gscat .* phi_inc;

    if calc_mua
      Gscat_mua = (exp(j * k * rd) ./ (4 * pi * rd) ...
		   - exp(j * k * rd1i) ./ (4 * pi * rd1i) ...
		   ) * volVoxel * pmiModel.v(idxLambda) / D;
    
      A(1:nPts, iMeas) = Gscat_mua .* phi_inc_mua;
    end    

    if calc_musp
      Gscat_musp_x = ((j*k -1./rd).*(Xm-Det(1)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
	  - (j*k -1./rd1i).*(Xm-Det(1)).*exp(j * k * rd1i) ./ (4 * ...
						  pi * rd1i.^2) ...
		      ) * volVoxel;
      
      Gscat_musp_y = ((j*k -1./rd).*(Ym-Det(2)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
	  - (j*k -1./rd1i).*(Ym-Det(2)).*exp(j * k * rd1i) ./ (4 * pi * rd1i.^2) ...
		      ) * volVoxel;
      
      Gscat_musp_z = ((j*k -1./rd).*(Zm-Det(3)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
		      - (j*k -1./rd1i).*(z1i-Det(3)).*exp(j * k * rd1i) ./ (4 * pi * rd1i.^2) ...
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
    ImageSrc = [Src(1) Src(2) 2*zBnd-Src(3)];
    rsrcdet = norm(Src - Det);
    rimgdet = norm(ImageSrc - Det);
    Phi_Inc(iMeas) = exp(j * k * rsrcdet) ./ (-4 * pi * rsrcdet) - ...
            exp(j * k * rimgdet) ./ (-4 * pi * rimgdet) ;

end

A = A.';
if Debug
    fprintf('\n');
end
