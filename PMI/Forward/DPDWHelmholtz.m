%DPDWHelmholtz  Forward DPDW model using exact solution of helmholtz equation.
%
%   [Phi_Inc] = DPDWHelmholtz(Model, MeasList, Debug)
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
%
%   Calls: Hlm2ptZB, Hlm2ptNB, Hlm2ptslab
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
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: DPDWHelmholtz.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.3  2000/04/13 20:23:17  dboas
%
%  Removed the display of 'Using david's defn of D'
%
%  Revision 3.2  1999/12/06 22:47:08  dboas
%  call to calcExtBnd needed correct lambda_index to the index of
%  refraction
%
%  Revision 3.1  1999/11/16 22:45:26  dab
%  Initial revision.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phi_Inc] = DPDWHelmholtz(Model, MeasList, Debug)

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
%%  Get the extrapolated boundary distance for cases of boundaries
%%
zBnd = calcExtBnd(Model.idxRefr(idxLambda), Model.Mu_sp(idxLambda));


if Debug
    fprintf('Modulation Freq: %f MHz\n', Model.ModFreq(idxFreq));
    fprintf('D = %e\n', D);
    fprintf('Re{k} = %f cm^-1\n', real(k));
    fprintf('Im{k} = %f cm^-1\n', imag(k));
    fprintf('Extrapolated boundary = %f cm\n', zBnd);
end



switch lower(Model.Boundary.Geometry)
  case { 'semi-infinite', 'semi', 'extrapolated'}
    if Debug
       fprintf(['Executing extrapolated zero boundary' ...
                       ' computation\n']);
    end
    [Phi_Inc] = Hlm2ptZB(Model, k, MeasList, zBnd, Debug);
          
  case {'infinite', 'inf'}
    if Debug
       fprintf(['Executing infinite medium boundary' ...
                       ' computation\n']);
    end
    [Phi_Inc] = Hlm2ptNB(Model, k, MeasList, Debug);
          
  case {'slab'}
    [Phi_Inc] = Hlm2ptSlab(Model, k, MeasList, zBnd, Debug);
        
  otherwise
    error(['Unknown boundary condition: ' Model.Boundary.Geometry]);
end

%%
%%  Scale A and Phi_Inc for the effective source amplitude in Helmholtz
%%  expression for the DPDW equation.  Addditionally, scale A as the
%%  mapping from del mu_a to the scattered fluence.  The Born-1
%%  approximation maps a change in k to the scattered fluence.
%%
Phi_Inc = -Model.v(idxLambda) / D * Phi_Inc;
