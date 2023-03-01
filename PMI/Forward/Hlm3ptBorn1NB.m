%HLM3PTBORN1NB  1st Born approx. of a perturbed Helmholtz equation, no boundary.
%
%   [A Phi_Inc] = Hlm3ptBorn1NB(pmiModel, k, MeasList, Debug)
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
%   Debug       OPTIONAL: Print out debugging info.
%
%   HLM3PTBORN1NB computes the forward weighting matrix associated with
%   the Born-1 approximation to a spatial varying k.  The incident response
%   at each voxel is first computed and then contribution each detector from
%   each voxel is computed.  This is done for each source as well.  The
%   sources are assumed to be unit amplitude point sources.
%
%   Do not specify the sources or detector exactly at a sampling point to
%   prevent dived by zero errors when evaluating the Green's functions.
%
%   Calls: none.
%
%   Bugs: none known.

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
%  $Date: 2001/02/02 17:52:37 $
%
%  $Revision: 1.4 $
%
%  $Log: Hlm3ptBorn1NB.m,v $
%  Revision 1.4  2001/02/02 17:52:37  dboas
%  Removed the 1/3 factor from the calculation of the perturbation matrix for
%  delta_D.  I don't recall
%  why this 1/3 appeared originally, but it shouldn't be there.
%
%  Revision 1.3  2000/11/07 16:46:56  dboas
%  Added v/D scaling for mua perturbation.
%  Added scattering perturbation.
%
%  Revision 1.2  2000/07/27 15:06:21  dboas
%  Added calculation for scattering perturbation.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.1  1999/11/15 16:22:41  rjg
%  Input parameter list changed to take pmiModel structure and MeasList table.
%
%  Revision 3.0  1999/06/17 17:39:56  rjg
%  Initial Revision for PMI 3.0
%
%  Revision 2.1  1999/02/05 20:55:07  rjg
%  Clearified help section and corrected the name of the function.
%
%  Revision 2.0  1998/08/05 15:50:55  rjg
%  Start of version 2
%  Handles structure for computational volume.
%
%  Revision 1.3  1998/06/18 15:00:47  rjg
%  Transformed x,y,z locations into column vectors so that it does not
%  have to be done later in functions of these variables.
%  Changed the variable Phi_Det to Phi_Inc since it is the incident field.
%
%  Revision 1.2  1998/06/05 18:12:56  rjg
%  Changed psi to phi.
%
%  Revision 1.1  1998/06/03 16:04:48  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = Hlm3ptBorn1NB(pmiModel, k, MeasList, Debug)

if nargin < 4
    Debug = 0;
end

nMeas = size(MeasList,1);
idxLambda = MeasList(1,4);

%%
%%  Extract the source and detector positions
%%
pSrc = getOptodePos(pmiModel.Src);
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
        fprintf('%d  ', iSrc);
    end
    Src = pSrc(MeasList(iMeas, 1),:);
    Det = pDet(MeasList(iMeas, 2),:);

    %%
    %%  Compute the incident response on every voxel from each measurement 
    %%  (source detector pair).  The source is represented by a point source.
    %%
    rs = sqrt((Xm - Src(1)).^2 + (Ym - Src(2)).^2 + (Zm - Src(3)).^2);
    
    if calc_mua
      phi_inc_mua = exp(j * k * rs) ./ (-4 * pi * rs);
    end
    if calc_musp
      phi_inc_musp_x = (j*k -1./rs).*(Xm-Src(1)).*exp(j * k * rs) ./ (-4 * pi * rs.^2);
    
      phi_inc_musp_y = (j*k -1./rs).*(Ym-Src(2)).*exp(j * k * rs) ./ (-4 * pi * rs.^2);
    
      phi_inc_musp_z = (j*k -1./rs).*(Zm-Src(3)).*exp(j * k * rs) ./ (-4 * pi * rs.^2);
    end

    %%
    %%  Compute the scattered response at this detector from all
    %%  of the voxels.  Note the that negative sign in the Green's function
    %%  is cancelled by the negative sign in the Born approximation
    %%
    rd = sqrt((Xm - Det(1)).^ 2 + (Ym - Det(2)).^2 + (Zm - Det(3)).^2);

    if calc_mua
      Gscat_mua = (exp(j * k * rd) ./ (4 * pi * rd) ...
		   ) * volVoxel * pmiModel.v(idxLambda) / D;
    
      A(1:nPts, iMeas) = Gscat_mua .* phi_inc_mua;
    end    

    if calc_musp
      Gscat_musp_x = ((j*k -1./rd).*(Xm-Det(1)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
		      ) * volVoxel;
      
      Gscat_musp_y = ((j*k -1./rd).*(Ym-Det(2)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
		      ) * volVoxel;
      
      Gscat_musp_z = ((j*k -1./rd).*(Zm-Det(3)).*exp(j * k * rd) ./ (4 * pi * rd.^2) ...
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
    rsrcdet = norm(Src - Det);
    Phi_Inc(iMeas) = exp(j * k * rsrcdet) ./ (-4 * pi * rsrcdet);

end
A = A.';
if Debug
    fprintf('\n');
end
