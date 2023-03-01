% addShotNoise  Add Poisson-distributed shot noise to simulated data
%
% Phi = addShotNoise(SD, MeasList, Phi0, gain, bwidth);
%
% Inputs:
%   Phi0 - previously generated synthetic data (in Watts)
%   gain - electronic gain of the imaging system (linear, not dB)
% bwidth - bandwidth of detector (in Hz)
%
% Outputs:
%   Phi  - Phi0 + noise
%

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

function [Phi] = addShotNoise(SD, MeasList, Phi0, gain, bwidth);

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

Phi = zeros(size(Phi0));

% Shot noise is wavelength-dependent, do each \lambda separately

for iWvl = 1:length(SD.Lambda)
  ml = find(MeasList(:,4) == iWvl);

  if (~isempty(ml))
     if (SD.Lambda(SD.MeasList(ml(1),4)) <= 0)
	error('Wavelength must be a positive number');
     end
     
     Phi(ml) = Poisson(Phi0(ml), gain, bwidth, SD.Lambda(iWvl) * 1e-9);
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Poisson distributed shot noise.  Phi is in Watts.
%
% I'm not convinced that real-imaginary is the right way to handle
% complex-valued data, but the alternative, amplitude-phase, is clearly
% wrong.
%
% Returns a vector Poisson-distributed random numbers with mean value
%  equal to the number of photons corresponding to the signal Phi0.

function [PhiN] = Poisson(Phi0, gain, bwidth, lambda)

if (lambda > 1e-6)
   % Assume looking at nm, but print warning message
   warning('Converting wavelength from nm to cm');
   lambda = lambda * 1e-9 * 100;
end

% Planck's constant times speed of light (J cm)
hc = 6.62606876e-34 * 2.99792458e10;

% Convert from power to number of photons
Nph = bwidth * abs(Phi0 / gain) / (hc / lambda);

% Generate Poisson-distributed random numbers and convert the deviation
% from the mean from photons back to powers.

if (isreal(Phi0))
   PhiN = poissrnd(max(0,Nph)) - max(0, Nph);
else
   Nph = Nph / sqrt(2);          % Split between I & Q
   
   PhiI = poissrnd(max(0,Nph)) - max(0, Nph);
   PhiQ = poissrnd(max(0,Nph)) - max(0, Nph);
   
   % Recombine with same noise on I and Q
   PhiN = PhiI + i*PhiQ;
   
   % Rotate to have the correct phase
   PhiN = PhiN .* exp(i*angle(Phi0)) * exp(-i*pi/4);
end

PhiN = gain * (PhiN * (hc / lambda) / bwidth);

% Each element of PhiN is a poisson-distributed random variable
%  with a mean of Phi0(k).

PhiN = PhiN + Phi0;

return;
