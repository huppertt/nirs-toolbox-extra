%HLM2PTZB  Exact solution of Helmholtz equation, z boundary.
%
%   [Phi_Inc] = Hlm2ptZB(pmiModel, k, MeasList, zBnd, Debug)
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
%   HLM2PTZB computes the exact solution of the helmholtz
%   equation for a spatial uniform k.  A zero fluence boudary
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
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: Hlm2ptZB.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.1  1999/11/16 22:53:09  dab
%  Initial revision.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phi_Inc] = Hlm2ptZB(pmiModel, k, MeasList, zBnd, Debug)

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

Phi_Inc = zeros(nMeas, 1) + j * ones(nMeas, 1);

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
    %%  Compute the incident response at the detector.
    %%
    ImageSrc = [Src(1) Src(2) 2*zBnd-Src(3)];
    rsrcdet = norm(Src - Det); 
    rimgdet = norm(ImageSrc - Det);
    Phi_Inc(iMeas) = exp(j * k * rsrcdet) ./ (-4 * pi * rsrcdet) - ...
            exp(j * k * rimgdet) ./ (-4 * pi * rimgdet) ;

end

if Debug
    fprintf('\n');
end
