%DPDWBorn1NB    Forward DPDW model using 1st born for a free space medium.
%
%   [A Phi_Inc] = DPDWBorn1NB(Model, MeasList, Debug)
%
%   A           The forward matrix relating the contribution from each voxel to
%               a specific source-detector pair.  Each row is for a different
%               source detector pair.
%
%   Phi_Inc     The incident response at each detector from each source in a
%               column vector.  The same combination pattern as the columns of A.
%
%   Model       The PMI Model structure contain the following fields: CompVol,
%               Mu_sp, Mu_a, v, idxRefr, and f.
%
%   MeasList    The measurement list for this simulation.  All frequencies and
%               wavelengths must be the same.
%
%   Debug       OPTIONAL: Print out debugging info.
%
%   DPDWBorn1B generate a Born-1 approximate forward model and the incident
%   field response for a single frequency and wavelength in the scenario
%   specified in Model.  MeasList specifies which sources and detectors are to
%   be used, only a single frequency and wavelength can be calculated per call
%   to this routine.
%   
%   Calls: Hlm3ptBorn1NB
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
%  $Log: DPDWBorn1NB.m,v $
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
%  Revision 1.3  2000/09/01 20:11:16  dboas
%  Moved calculation of coupling coefficient matrix to genFwdMat.m
%
%  Added Model.Method.Type 'BornN' and 'FullBorn'.
%
%  Revision 1.2  2000/08/01 20:06:07  dboas
%  Added calculation of matrix to solve for S D coupling coefficients.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.2  1999/11/09 22:24:57  rjg
%  Changed calling format.  It now handles the processing in MeasList as long
%  as MeasList is a single frequency and wavelength (constant k).
%
%  Revision 3.1  1999/09/29 21:31:45  rjg
%  Fixed case problem for function calls under unix.
%
%  Revision 3.0  1999/06/17 17:36:01  rjg
%  Initial version for PMI 3.0
%
%  Revision 2.1  1999/02/05 20:53:04  rjg
%  Again, changed the value of D, it is now
%  v / (3 * (mu_sp + mu_a))
%
%  Revision 2.0  1998/09/11 21:27:08  rjg
%  Start of version 2
%  Handles structure for computational volume
%
%  Revision 1.4  1998/06/18 15:07:31  rjg
%  Changed D again to match David's D.  1. removed mu_a from D and 
%  2. put v back into D, this did not change k since it was accounted for
%  in previous versions.
%  Scaling of Phi_Inc and A from the Born-Helmholtz approximation was
%  modified for a D without mu_a and negative sign in A that was dropped has
%  be accounted for.
%
%  Revision 1.3  1998/06/10 18:33:27  rjg
%  D is now in units of length, v was removed and move into the imaginary
%  part of k.  This is the same as dividing through the diffusion equation
%  by v.
%  Both Phi_Inc and A are now scaled correctly for the Born-1 approximation,
%  see derivation in notes.
%
%  Revision 1.2  1998/06/05 18:13:39  rjg
%  Changed psi to phi.
%
%  Revision 1.1  1998/06/03 16:03:16  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = DPDWBorn1NB(Model, MeasList, Debug)

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

if Debug
    fprintf('Modulation Freq: %f MHz\n', Model.ModFreq(idxFreq));
    fprintf('D = %e\n', D);
    fprintf('Re{k} = %f cm^-1\n', real(k));
    fprintf('Im{k} = %f cm^-1\n', imag(k));
end

%%
%%  Calculate the Born approximation and incident fields
%%
if strcmpi(Model.Method.Type,'Born') | strcmpi(Model.Method.Type,'Rytov')
  [A Phi_Inc] = Hlm3ptBorn1NB(Model, k, MeasList, Debug);
  Phi_Inc = -Model.v(idxLambda) / D * Phi_Inc;
%  A = (Model.v(idxLambda) / D)^2 * A;
  A = (Model.v(idxLambda) / D) * A;
elseif strcmpi(Model.Method.Type,'FullBorn') | ...
      strcmpi(Model.Method.Type,'BornN') | ...
      strcmpi(Model.Method.Type,'ExtBorn') | ...
      strcmpi(Model.Method.Type,'ExtBornN')
  [A Phi_Inc] = HlmFullBorn_NB(Model, k, MeasList, Debug);
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
	
	
