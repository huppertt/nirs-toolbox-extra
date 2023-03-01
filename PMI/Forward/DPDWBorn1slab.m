%DPDWBorn1slab    Forward DPDW model using 1st born for a slab medium.
%
%   [A Phi_Inc] = DPDWBorn1slab(Model, MeasList, Debug)
%
%   A           The forward matrix relating the contribution from each voxel to
%               a specific source-detector pair.  Each row is for a different
%               source detector pair.
%
%   Phi_Inc     The incident response at each detector from each source
%               in a column vector.  The same combination pattern as the
%               columns of A.
%
%   Model       The PMI Model structure contain the following fields: CompVol,
%               Mu_sp, Mu_a, v, idxRefr, and f.
%
%   MeasList    The measurement list for this simulation.  All frequencies and
%               wavelengths must be the same.
%
%   Debug       OPTIONAL: Print out debugging info.
%
%   The true z boundary is assumed to be at z=0.
%
%   Calls: Hlm3ptBorn1slab
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
%  $Date: 2001/02/07 16:30:15 $
%
%  $Revision: 1.6 $
%
%  $Log: DPDWBorn1slab.m,v $
%  Revision 1.6  2001/02/07 16:30:15  dboas
%  Added forward model method type 'ExtBornN'.  This does an iterative extended
%  Born approximation.
%  DON'T FORGET THAT 'BornN' AND 'ExtBornN' DON'T NECESSARILY CONVERGE IF THE
%  PERTURBATION IS TOO LARGE.
%
%  Revision 1.5  2001/02/02 17:52:13  dboas
%  Added Extended Born method
%
%  Revision 1.4  2000/11/07 16:46:38  dboas
%  The matrix A was scaled by v/D specific for mua perturbations.  This was moved
%  to a lower level because of the potential for scattering perturbations which
%  are not scaled by v/D.
%
%  Revision 1.3  2000/09/01 20:13:11  dboas
%  Moved calculation of coupling coefficient matrix to genFwdMat.m
%
%  Added Model.Method.Type 'BornN' and 'FullBorn'.
%
%  Revision 1.2  2000/07/27 15:05:32  dboas
%  Scale the A matrix by (v/D).^2 .
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1999/12/06 22:46:11  dboas
%  call to calcExtBnd needed correct lambda_index to the index of
%  refraction
%
%  Revision 1.1  1999/11/18 14:18:55  tgaudett
%  Creation of Slab Geometry Calculation
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = DPDWBorn1slab(Model, MeasList, Debug)

if nargin < 3
    Debug = 0;
end
%%
%%  Get the frequency and wavelength index from MeasList
%%
idxFreq = MeasList(1,3);
if ~all(idxFreq == MeasList(:,3))
    error('all frequency indices must be the same in MeasList')
end
idxLambda = MeasList(1,4);
if ~all(idxLambda == MeasList(:,4))
    error('all wavelength indices must be the same in MeasList')
end

%%
%%  Compute the PDE parameters from the media parameters
%%
%D = v / (3 * (mu_sp + mu_a));
%disp('Using David''s defn: D = v / (3 * mu_sp)')
D = Model.v(idxLambda) / (3 * Model.Mu_sp(idxLambda));
k = sqrt(-Model.v(idxLambda) * Model.Mu_a(idxLambda) / D + ...
    j * 2 * pi * Model.ModFreq(idxFreq) * 1E6 / D);

%%
%%  Get the extrapolated boundary distance
%%
zBnd = calcExtBnd(Model.idxRefr(idxLambda), Model.Mu_sp(idxLambda));

if Debug
    fprintf('Modulation Freq: %f MHz\n', Model.ModFreq(idxFreq));
    fprintf('D = %e\n', D);
    fprintf('Re{k} = %f cm^-1\n', real(k));
    fprintf('Im{k} = %f cm^-1\n', imag(k));
    fprintf('Extrapolated boundary = %f cm\n', zBnd);
end

%%
%%  Calculate the Born approximation and incident fields
%%
if strcmpi(Model.Method.Type,'Born') | strcmpi(Model.Method.Type,'Rytov')
  [A Phi_Inc] = Hlm3ptBorn1slab(Model, k, MeasList, zBnd, Debug);
  Phi_Inc = -Model.v(idxLambda) / D * Phi_Inc;
%  A = (Model.v(idxLambda) / D)^2 * A;
  A = (Model.v(idxLambda) / D) * A;
elseif strcmpi(Model.Method.Type,'FullBorn') | ...
      strcmpi(Model.Method.Type,'BornN') | ...
      strcmpi(Model.Method.Type,'ExtBorn') | ...
      strcmpi(Model.Method.Type,'ExtBornN')
  [A Phi_Inc] = HlmFullBorn_slab(Model, k, MeasList, zBnd, Debug);
  Phi_Inc = -Model.v(idxLambda) / D * Phi_Inc;
%  A = -(Model.v(idxLambda) / D) * A;
  A = - A;
end

%%
%%  If Method = 'Rytov' then normalize the matrix by PhiInc
%%
if strcmp(Model.Method.Type,'Rytov')
  A = A ./ (Phi_Inc * ones(1,size(A,2)));
end
	


      
      