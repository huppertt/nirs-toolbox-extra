% FD3PT   Generate the incident field for FD model (1st Born)
%
% [Phi0,A] = FD3pt(SD, Medium, MeasList, OptProp, [Debug]);
%   
%   SD, Medium - PMI data structures
%   MeasList   - The measurement list
%   OptProp    - 2x1 vector of matrix flags.
%                OptProp(1) -> allow absorbing perturbations
%                OptProp(2) -> allow scattering perturbations
%
% Returns:
%    Phi0 - The total measured fluence for each
%         source-detector pair specified in SD.MeasList given a
%         homogneeous medium.  
%    A    - The forward matrix.

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

function[Phi0,A] = FD3pt(SD, Medium, MeasList, OptProp, Debug)

if (~exist('MeasList','var')) | (isempty(MeasList))
   MeasList = SD.MeasList;
end

if ~exist('Debug','var')
   Debug = 0;
end

switch lower(Medium.Geometry)
   case {'infinite', 'inf', 'inft' }
      if Debug
	 fprintf(['Executing infinite medium boundary computation\n']);
      end
      
      [Phi0,A] = Hlm3ptBorn1NB(SD, Medium, MeasList, OptProp, Debug);
      
   case { 'semi-infinite', 'semi', 'extrapolated'}
      if Debug
	 fprintf(['Executing extrapolated zero boundary computation\n']);
      end

      [Phi0,A] = Hlm3ptBorn1ZB(SD, Medium, MeasList, OptProp, Debug);
      
   case {'slab'}
      if Debug
	 fprintf(['Executing slab boundary computation\n']);
      end

      [Phi0,A] = Hlm3ptBorn1SB(SD, Medium, MeasList, OptProp, Debug);
      
   otherwise
      error(['Unknown or unsupported boundary condition: ' Medium.Geometry]);
end

% The default NMeas by NVox packing is wrong.  The matrix should have a
% block structure so that optical properties at different wavelengths can't
% talk to each other.  Fix this here for now, go into the individual
% routines to fix it when I get a chance.

wl = unique(MeasList(:,4));
nw = length(SD.Lambda);

if (nw > 1)
   % If there's only one wavelength, the measurement happens to have the
   % correct packing.  Repack if there's more than one wavelength.
   
   [nm,nv] = size(A);
   
   % Matlab does a miserable job copying from full into sparse matrices.
   
   if (exist('repackBornFwdMat') == 3)
      % The fast way, assumes compiler support
      
      Ap = repackBornFwdMat(SD, MeasList, A, OptProp, Debug);
      
      if (Debug)
	 disp('Finished repacking');
      end
   else
      % The slow way, will always work
      
      % I happen to know the exact size of the matrix, take advantage of
      % that here and call spalloc instead of sparse(zeros(..)).

      if (Debug)
	 disp('Block packing forward matrix');
      end

      if all(OptProp(1:2) ~= 0)
	 Ap = spalloc(nm, nw*(2*nv), nm*(2*nv));
      else
	 Ap = spalloc(nm,     nw*nv,     nm*nv);
      end
      
      for iwvl = 1:length(SD.Lambda)
	 if (Debug)
	    disp([ '  ...Wavelength ' num2str(iwvl) ]);
	 end
	    
	 llist = find(MeasList(:,4) == iwvl);
	 
	 if (~isempty(llist))
	    for ivox = 1:nv
	       ivoxp = (iwvl-1)*nv + ivox;
	       
	       Ap(llist,ivoxp) = A(llist,ivox);
	    end
	 end
      end
   end

   clear A;
   A = Ap;
   clear Ap;
end

return;
