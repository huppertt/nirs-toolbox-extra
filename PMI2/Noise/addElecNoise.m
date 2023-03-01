% addElecNoise  Add Gaussian-distributed electronic noise to simulated data
%
% Phi = addElecNoise(SD, MeasList, Phi0, SNR);
%
% Inputs:
%   SD, MeasList - not used, for uniformity with addShotNoise
%
%   Phi0 - previously generated synthetic data.  Can also be used
%          to specify a noise-equivalent power if the SNR is set to 1.
%   SNR  - average signal to noise ratio (linear, NOT dB!)
%
% Outputs:
%   Phi  = Phi0 + noise

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

function [Phi] = addElecNoise(SD, MeasList, Phi0, SNR);

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Phi0','var'))
   Phi0 = [];
end

if (~exist('SNR','var') | isempty(SNR))
   % Phi0 can also be used to specify the NEP
   SNR = 1;
end

if (length(Phi0) == 1)
   % Phi0 has been used to specify the NEP instead of the signal.  Use this
   % value to generat an instance of noise with the user-requested NEP.
   
   NEP = Phi0;
   
   Phi = NEP * randn(size(MeasList,1),1);
else
   Phi = Phi0 + Gaussian(SNR, Phi0);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the Gaussian distributed noise

function [N] = Gaussian(SNR, Phi0)

signal = sqrt(mean(abs(Phi0).^2));      % Define signal to be RMS signal
noise  = signal / SNR;

repack = 0;

if (~isreal(Phi0))
  [n1,n2] = size(Phi0);

  % Break Phi0 into real and imaginary components
  Phi0 = [ real(Phi0(:)), imag(Phi0(:)) ];
  repack = 1;
end

% Noise
N = noise * randn(size(Phi0));

if (repack)
   % Pack back to complex numbers

   N = N(:,1) + i*N(:,2);
   N = reshape(N(:), n1, n2);
end

return;
