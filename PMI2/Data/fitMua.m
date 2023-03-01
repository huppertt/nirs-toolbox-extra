% fitMua   Fit optical absorption to data
%
% Mua = fitMua(SD, Medium, MeasList, data);
%
% Inputs:
%   SD, Medium - PMI2 structures
%   data       - Measured fluence vector.
%
% Outputs:
%   Mua        - Fit of Mua.  For wavelengths not in the measurement
%                list, Mua is copied from Medium.Muao.

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

function[Mua] = fitMua(SD, Medium, MeasList, Phi, Debug);

% Copy out of SD if missing or not provided

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

% Default values
Mua = Medium.Muao;

slst  = unique(MeasList(:,4));
hcost = @muaCost;

opts  = optimset('fminbnd');

for k = 1:length(slst)
   iWvl = slst(k);

   ml = find(MeasList(:,4) == iWvl);

   Mua(iWvl) = fminbnd(hcost, 0.005, 0.200, opts, ...
		       SD, Medium, MeasList(ml,:), Phi(ml), Debug);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cost = tmpcost(mu, SD, Medium, MeasList, Phi)

function[cost] = muaCost(mu, SD, Medium, MeasList, Phi, Debug)

slst = unique(MeasList(:,4));

Medium.Muao(slst) = mu;

if (isempty(Phi))
   warning('Empty data vector');
   keyboad;
end

% Get incident fluence, don't need full forward matrix for chisqr

Phi0 = DPDWHelmholtz(SD, Medium, MeasList);

if (isreal(Phi0) & (min(Phi0) < 0))
   warning('genBornData() returned a negative fluence');
   keyboard;
end

% Rescale the (possibly corrected) incident fluence

[sA, dA, pA] = fitSD(SD, Medium, MeasList, Phi);
Phi0         = Phi0 .* pA;

if (Debug)
   dR = calcSep(SD, MeasList);
   [dr,Idx] = sort(dR);
   
   figure(1);
   semilogy(dr, abs(Phi0(Idx) ./ pA(Idx)), dR, abs(Phi ./ pA),'ro')
   drawnow;
end

if (0)
   % Compute chi-squared for this set of parameters, for Poisson the
   % variance of the signal is directly proportional to the signal

   cost = sum(abs(Phi - Phi0).^2 ./ abs(Phi)) / (length(Phi)-1);
else
   % Use a Rytov-like cost function that weights each point
   %  more-or-less evenly
   
   cost = sqrt(mean(log(abs(Phi./Phi0)).^2));
end

if (Debug)
   disp(sprintf('  ..chisqr(%f) = %e', mu, cost));
end

return;
