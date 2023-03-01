% PMI_Example - sample PMI code, basic introduction to the toolbox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2004, Jonathan Stott
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tell matlab where to find the PMI files.  This should be done only
% once per session; I recommend putting it into your startup.m file.

% pmipath('/homes/monte/1/home/jstott/matlab/PMI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a few variables for convenience
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doscat = 1;		% flags, 1 or 0
doabs  = 1;
thick  = 6;		% Slab thickness

Method = 'Rytov';       % Born or Rytov

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate SD structure with source/detector/measurement information

SD.Lambda  = [690 800]; % Define 2 Wavelengths
SD.ModFreq = [  100  ]; % Single-frequency RF imager

% Source and Detector Positions

SD.SrcPos = SetOptode([-2:2:2], [-2:2:2], 0, 1);
SD.DetPos = SetOptode([-3:2:3], [-3:2:3], thick, 1);

% Source and Detector Amplitudes

SD.SrcAmp = 1e-3*ones(size(SD.SrcPos,1), ...
                      length(SD.Lambda), length(SD.ModFreq));
SD.DetAmp = 1e-3*ones(size(SD.DetPos,1), ...
                      length(SD.Lambda),length(SD.ModFreq));

% Generate a measurement list with all possible measurements

SD.MeasList = genMeasList(SD, 'all');

nWvl = length(SD.Lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate structure describing Medium to be imaged

% Optical properties

Medium.idxRefr = [ 1.37 1.36 ];
Medium.Muao    = [ 0.02 0.04 ];
Medium.Muspo   = [ 9.50 8.35 ];

% Imaging volume

Medium.Geometry       = 'slab';
Medium.Slab_Thickness = thick;

% Volume to reconstruct

Medium.CompVol.Type  = 'computed';
Medium.CompVol.XStep = 0.5;
Medium.CompVol.X     = [ -4.0 4.0 ];
Medium.CompVol.YStep = 0.5;
Medium.CompVol.Y     = [ -4.0 4.0 ];
Medium.CompVol.ZStep = 0.5;	% Avoid the problematic z=0 case
Medium.CompVol.Z     = [ 0.25 thick ];

nVox = length(sampleVolume(Medium.CompVol));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate forward matrix ('Born'/'Rytov' solve the homogeneous problem)

% Two wavelengths were defined, but I'm only going to reconstruct data
%  from one of them, as an example of using a subset of the measurement
%  list (and to save time).  Also drop the long-distance imaging pairs.

dR = calcSep(SD, SD.MeasList);

ml1 = find(SD.MeasList(:,4) == 1 & dR < 9);
MeasList = SD.MeasList(ml1, :);
clear ml1 dR;

disp('Generating forward matrix - please wait');

[Phi0, A] = genBornMat(SD, Medium, MeasList, Method, [doabs doscat]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define perturbations - NOTE: abs/scat perturbations supported if and
%                         only if they're supported by the forward
%                         problem (i.e., they're included in A)

if (doabs)
   Medium.Object{1}.Type   = 'Sphere';         % Spherical scatterer
   Medium.Object{1}.Pos    = [ 1 1.25 2.9 ];
   Medium.Object{1}.Radius = 1.2;
   Medium.Object{1}.Mua    = Medium.Muao;
   Medium.Object{1}.Musp   = Medium.Muspo - 1;
end

if (doscat)
   if (doabs)
      n = 2;
   else
      n = 1;
   end
   
   Medium.Object{n}.Type   = 'Sphere';          % Spherical absorber
   Medium.Object{n}.Pos    = [ -2.25 -1.25 1.75 ];
   Medium.Object{n}.Radius = 1.2;
   Medium.Object{n}.Mua    = Medium.Muao + 0.02;
   Medium.Object{n}.Musp   = Medium.Muspo;
   
   clear n;
end

disp('Generating simulated data');

Phi = genBornData(SD, Medium, MeasList, Method, Phi0, A, [doabs doscat]);

disp('Adding noise');

Phi = addShotNoise(SD,MeasList,Phi,10,1);    	% approx 1000:1 SNR
Phi = addElecNoise(SD,MeasList,Phi,1000);    	% average 1000:1 SNR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Invert the simulated data

disp('Inverting simulated data');

% A gives me d\mu_a and dD/D, convert to d\mu_a/\mu_a and dD/D to make
% both halves of order 1

if (doabs & doscat)
   for k = 1:nWvl
      % First half at each wavelength is absorbtion, second scattering
      
      A(:,2*(k-1)*nVox + [1:nVox]) = ...
	  A(:,2*(k-1)*nVox + [1:nVox]) * Medium.Muao(k);
   end
end

% Repack as real and imaginary for inversion

A1 = [ real(A); imag(A) ];

if (strcmpi(Method,'Rytov'))
   Y = [ real(log(Phi./Phi0)); imag(log(Phi./Phi0)) ];
else
   Y = [ real(Phi-Phi0); imag(Phi-Phi0) ];
end

clear A Phi0;

% Invert the data using filtered back-projection with Tikhonov
% regularization.  The regularization parameter alpha was adjusted 
% manually to get good images.

alpha = 1e-6 * normest(A1).^2;

X = fbp(A1, Y, alpha);

if (doabs & doscat)
   X = reshape(X, nVox, 2*nWvl);

   for k = 1:nWvl
      X(:,2*k-1) = X(:,2*k-1) *  Medium.Muao(k);
      X(:, 2*k ) = X(:, 2*k ) * -Medium.Muspo(k) / 3;
   end
else
   X = reshape(X, nVox, nWvl);
   
   if (doscat)
      for k = 1:nWvl
	 X(:,k) = X(:,k) * -Medium.Muspo(k) / 3;
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the reconstructed data [d\mu_a and d\mu_s']

if (doabs & doscat)
   dmua = X(:, 1:2:end);
   dmus = X(:, 2:2:end);
elseif (doabs)
   dmua = X;
   dmus = zeros(0,nWvl);
else
   dmua = zeros(0,nWvl);
   dmus = X;
end

clear X;

for k = 1:nWvl
   if (any(MeasList(:,4) == k))
      showImage(Medium, dmua(:,k), dmus(:,k));
   end
end

% Display the actual scattering perturbation(s), just as a reference

for k = 1:nWvl
   % Only display wavelengths actually used
   
   if (any(MeasList(:,4) == k))
      dmua0 =  calcDelMuA(Medium, [], k);
      dmus0 = calcDelMuSp(Medium, [], k);

      if (~isempty(dmua0) & all(dmua0(:) == 0))
	 dmua0 = [];
      end

      if (~isempty(dmus0) & all(dmus0(:) == 0))
	 dmus0 = [];
      end

      showImage(Medium, dmua0(:), dmus0(:));
   end
end

% Clear some variables before exiting to save memor

clear dmu0 k n nWvl nVox thick k alpha;
