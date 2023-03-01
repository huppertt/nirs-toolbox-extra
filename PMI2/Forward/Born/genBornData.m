% GENBORNDATA  Generate noiseless measured data using the specified method
%
% [Phi, Phi2pt] = genBornData(SD, Medium, MeasList, Method
%                             Phi0, A, OptProp, MethodParam1 );
%
%   Method - This specifies the method used to solve the Forward
%            Problem.  It can be one of the following:
%
%               'Helmholtz Homogeneous' - Solves the analytic
%                  Green's function for the geometry specified in
%                  Medium.Geometry.
%
%               'Born' - Solves the first Born approximation.
%
%               'Rytov' - Solves the Rytov approximation.
%
%               'BornN' - Solve the Nth Born approximation.
%                  MethodParam1 - The number of iterations N
%                                 to perform
%
%   SD, Medium - The PMI data structures
%   MeasList   - The Measurement List.  For "FullBorn" method, the 
%                measurement list really should match the measurement 
%		 list used to generate A and Phi0.
%
%   Phi0 - The incident fluence, not required for the method
%          'Helmholtz Homogeneous'.
%
%   A    - The forward matrix, not required for the Method 
%          'Helmholtz Homogeneous'.
%
%   OptProp - is a flag with 2 elements. 1-Yes, 0-No.
%      Default = [1 1]  (OptProp(1)->abs, OptProp(2)->scat)
%
%   MethodParam1 - A parameter which is Method dependent (n for nth-Born,
%        radius for extended born). See code for details.
%
% Returns
%   Phi    - the total measured fluence
%   Phi2pt - the total fluence at the voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2004, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang,
%                     Jonathan Stott
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

function  [Phi, Phi2pt] = genBornData(SD, Medium, MeasList, Method, ...
				      Phi0, A, OptProp, MethodParam1)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

% Initialize to empty lists

Phi    = [];
Phi2pt = [];

% No objects, no work

if (~isfield(Medium,'Object'))
   warning('Medium.Object not defined, total field equals incident');
   keyboard;
   
   Phi = Phi0;
end

%%
%%  Generate measured data
%%

if (length(OptProp) < 2)
   OptProp(end:2) = 0;
end

% Get perturbation flags, make sure either
%  1) calc_mua or calc_musp (or both) are non-zero
%  2) calling method is "Helmholtz Homogeneous" (i.e. unperturbed)

calc_mua = OptProp(1);
calc_musp= OptProp(2);

if (strcmpi(Method, 'Helmholtz Homogeneous') ~= 0)
   if ((calc_mua == 0) & (calc_musp == 0))
      error('OptProp() cannot be all zeros');
   end
end

switch Method
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% HELMHOLTZ HOMOGENEOUS (analytical solution, unperturbed)
   %%
   case 'Helmholtz Homogeneous'
      
      Phi0 = DPDWHelmholtz(SD, Medium, MeasList);
	 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% BORN and RYTOV approximations
   %%
   case { 'Born', 'Rytov' }
      %%
      %%  Generate the volume data using the 1st Born method(s)
      %%

      delMu   = zeros(size(A,2),1);

      for iWvl = 1:length(SD.Lambda)
	 lMeas   = find(MeasList(:,4) == iWvl);

	 if (~isempty(lMeas))
	    [dMu, dD] = getMuPert(Medium, iWvl, calc_mua, calc_musp);

	    % Voxels change faster than wavelengths
	 
	    if ((calc_mua & isempty(dMu)) | (calc_musp & isempty(dD)))
	       error('Perturbation requested, but result is empty vector');
	    end
	    
	    if (~isempty(dMu) & ~isempty(dD))
	       nvox = length(dMu);
	    
	       delMu((iWvl-1)*(2*nvox)        + [1:nvox]) = dMu(:);
	       delMu((iWvl-1)*(2*nvox) + nvox + [1:nvox]) = dD(:);
	    elseif (~isempty(dMu))
	       nvox = length(dMu);
	       delMu((iWvl-1)*nvox + [1:nvox]) = dMu(:);
	    elseif (~isempty(dD))
	       nvox = length(dD);
	       delMu((iWvl-1)*nvox + [1:nvox]) = dD(:);
	    else
	       error('Impossible state: delMua and delD both empty');
	    end
	 end
      end

      PhiScat = A * delMu(:);

      % First Born and First Rytov approximations are very closely related
	 
      if (strcmpi(Method,'Born'))
	 Phi = Phi0 + PhiScat;
      elseif (strcmpi(Method, 'Rytov'))
	 Phi = Phi0 .* exp(PhiScat);
      else
	 error([ 'Unknown method ' Method ]);
      end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% EXTENDED BORN
   %%
   %% Similar to N-th Born and Full Born, but not yet implemented.
   %%
   case { 'ExtBornN', 'ExtBorn' }
      error('Extended Born not implemented in v2.0 yet');
   
      % Integrate code later
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% N-th BORN's, FULL BORN, and EXTENDED BORN
   %%
   %% Most of the calculation is the same, so lump those parts together
   %% and vary the calculation only at the last step.
   %% 
   %% Extended born not implemented yet.
   case { 'BornN', 'FullBorn' }

      % Recreate the mapping used to generate A and Phi0.  This will
      % likely fail if the measurement list has changed since A and
      % Phi0 were generated, but I don't see any way around that.

      FB = MLtoFB(SD, Medium, MeasList, OptProp);

      % Allocate space for measured fluence
      Phi = zeros(size(MeasList,1),1);
      
      % Really \delta D / D and dMu_a, not dMu_s^' and dMu_a
      delA = zeros(length(FB.optF) * FB.nVox, 1);
      delS = zeros(length(FB.optF) * FB.nVox, 1);  
      
      for iOpt = 1:length(FB.optF)
	 iwvl = MeasList(FB.optF(iOpt),4);
	 
	 ml = find(MeasList(:,4)==iwvl);
	 
	 if (isempty(ml))
	    % Nothing to do
	    continue;
	 end
	 
	 [delMua, delD] = getMuPert(Medium, iwvl, calc_mua, calc_musp);

	 if ((calc_mua & (length(delMua(:)) ~= FB.nVox)) | ...
	     (calc_musp& (length(delD(:))   ~= FB.nVox)))
	    warning('FullBorn and Perturbation have different voxel sizes');
	    keyboard;
	 end

	 if (calc_mua)
	    delA((iOpt-1)*FB.nVox + [1:FB.nVox]) = ...
                delA((iOpt-1)*FB.nVox + [1:FB.nVox]) + delMua(:);
	 end
	 
	 if (calc_musp)
	    delS((iOpt-1)*FB.nVox + [1:FB.nVox]) = ...
                delS((iOpt-1)*FB.nVox + [1:FB.nVox]) + delD(:);
	 end
      end

      % Pack down perturations to match forward matrix
      
      for k = 1:FB.nOpt
         v1 = (iOpt - 1)*FB.nVox + [1:FB.nVox];
         v2 = ( k   - 1)*FB.nVox + [1:FB.nVox];
         
         delA(v2) = delA(v1);
         delS(v2) = delS(v1);
      end
      
      delA = delA(1:FB.nOpt*FB.nVox);
      delS = delS(1:FB.nOpt*FB.nVox);
      
      % Now that I have the perturation at all wavelengths, use the
      % perturbations to scale the forward problem.
      
      if (issparse(A))
	 % The rowscale()'s _really_ don't like sparse matrices.  
         % Unfortunately, my system memory _really_ doesn't like 
         % full matrices either.
         
	 A = full(A);
      end

      if (calc_mua & calc_musp)
	 % Rescale the voxel->voxel terms
         nvox = 4*FB.nVox*FB.nOpt;

	 A(:,1:nvox) = rowscale(A(:,1:nvox), ...
                                [ delA(:); delS(:); delS(:); delS(:) ]);
      elseif (calc_mua)
	 % Rescale the voxel->voxel terms
         nvox = 1*FB.nVox*FB.nOpt;

	 A(:,1:nvox) = rowscale(A(:,1:nvox), delA);
      elseif (calc_musp)
	 % Rescale the voxel->voxel terms
         nvox = 3*FB.nVox*FB.nOpt;
         
	 A(:,1:nvox) = rowscale(A(:,1:nvox), [ delS(:); delS(:); delS(:) ]);
      end

      if (strcmpi(Method,'FullBorn'))
	 %%
	 %% FULL BORN
	 %%
	    
         if (calc_musp)
            % The problem is really one of inverting non-square matrices.
	    error('Full Born not implemented for scattering perturbations');
         end
	 
	 % Unfortunately, "\" is not implemented for sparse complex
         % matrices, so A needs to be full to compute PhiTot.
	 
         if (issparse(A))
	    A = full(A);
	 end
	 
	 % Based on relation: 1 + x + x^2 + ... = (1 - x)^{-1} for |x| < 1
	 
	 % Left-divide only works for square matrices, which means this is 
         % only valid for absorbtion-only perturbations.
         
	 PhiTot = (speye(size(A)) - A) \ Phi0;

      elseif (strcmpi(Method,'BornN'))
	 %%
	 %% N-th BORN
	 %%
	    
         if (calc_musp)
	    error('N-th Born not implemented for scattering perturbations');
         end
	 
	 opts.display = 0;
	 maxeig = abs(eigs(A, 1, 'LM', opts));
     
	 if (maxeig > 1)
	    disp(sprintf('*** WARNING: Largest eigenvalue %f > 1', maxeig));
	    warning('*** BornN matrix will not converge!');
	    keyboard
	 end
	    
	 Pert   = Phi0;
	 PhiTot = Phi0;
	    
	 for n = 1:MethodParam1
	    if (calc_musp)
 	       % Pert_Grad Handles both calc_mua=0 and calc_mua=1
	       
 	       Pert_old = Pert;
 	       Pert     = A * Pert_Grad(Pert, Medium, calc_mua);
	       
 	       if max(abs(Pert(:))) > max(abs(Pert_old(:)))
 		  disp(sprintf('*** WARNING: Largest eigenvalue %f > 1', ...
 			       maxeig));
 		  warning('*** BornN matrix will not converge!');
 		  keyboard
 	       end
 	    elseif (calc_mua)
	       Pert = A * Pert;
	    end
	    
	    PhiTot = PhiTot + Pert;
	 end
      else
	 error([ 'Unknown method ' Method ]);
      end

      %%
      %% CONVERT BACK TO MEASLIST ORDER
      %%
	 
      vSrc = FB.srcR;
      vDet = FB.detR;

      % Accessing single elements is faster for sparse matrices,
      % vector operations are faster for full matrices
	 
      if (issparse(PhiTot))
	 for k = 1:length(vSrc)
	    Phi(k) = PhiTot(FB.nVox*FB.nOpt + vDet(k), vSrc(k));
	 end
      else
	 Phi = diag(PhiTot(FB.nVox*FB.nOpt + vDet, vSrc));
      end

      % EXTRACT THE 2PT FLUENCE WITHIN THE MEDIUM

      Phi2pt = PhiTot(1:FB.nVox*FB.nOpt, :);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% ERROR, UNKNOWN METHOD
   %%
   otherwise
      error(['Unknown forward method: ' Method])
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Turn Medium.Object{}'s into two perturbation vectors (scattering and
% absorption).  Objects may be an empty list (in which case
% Medium.Objects is used directly).
%
% All the different methods use this code, so it made more sense to
% put it into its own function.

function[delMua, delD] = getMuPert(Medium, iWvl, calc_mua, calc_musp)

if (calc_mua)
   delMua = calcDelMuA(Medium, Medium.Object, iWvl);
   delMua = delMua(:);
else
   delMua = [];
end

% Really \delta D / D, despite the name of the variable

if (calc_musp)
   delMu_sp = calcDelMuSp(Medium, Medium.Object, iWvl);
   delD     = -3 & delMu_sp(:);
else
   delD     = [];
end

return;


