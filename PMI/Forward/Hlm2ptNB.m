%HLM2PTZB  Exact solution of Helmholtz equation, infinite medium.
%
%   [Phi_Inc] = Hlm2ptZB(pmiModel, k, MeasList, Debug)
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
%   HLM2PTNB computes the exact solution of the helmholtz
%   equation for a spatial uniform k for an infinite medium.
%   The sources are assumed to be unit
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
%  $Log: Hlm2ptNB.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.1  2000/04/20 14:07:04  dboas
%
%  Initial version
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Phi_Inc] = Hlm2ptNB(pmiModel, k, MeasList, Debug)

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
    rsrcdet = norm(Src - Det);
    Phi_Inc(iMeas) = exp(j * k * rsrcdet) ./ (-4 * pi * rsrcdet);

end

if Debug
    fprintf('\n');
end
