% Validate_TD  Script originally used to validate the time-domain.
%              Now serves just as more sample code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear SD Medium

%% Basic parameters

Thickness = 5;
dX        = 2.5;
dY        = 1.5;
dZ        = 3.0;
domua     = 1;
domus     = 1;

%% Build the Medium structure

Medium.Muspo   = 8;
Medium.Muao    = 0.05;
Medium.idxRefr = 1.4;
Medium.v       = 2.99e10./Medium.idxRefr;

%Medium.Geometry       = 'infinite';
%Medium.Geometry       = 'semi';
Medium.Geometry       = 'slab';

Medium.Slab_Thickness = Thickness;

Medium.CompVol.Type  = 'uniform';
Medium.CompVol.XStep = 0.5;
Medium.CompVol.X     = dX;
Medium.CompVol.YStep = 0.5;
Medium.CompVol.Y     = dY;
Medium.CompVol.ZStep = 0.5;
Medium.CompVol.Z     = dZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build a frequency-domain SD structure and solve the forward problem

disp('Calculating FD problem');

clear SD

SD.Lambda = 800;

SD.SrcPos = [ 0 0 0 ];
SD.DetPos = [ 4 0 0 ];
SD.SrcAmp(1,1) = 1;
SD.DetAmp(1,1) = 1;
SD.ModFreq = [0:1023 -1024:-1]/(2048*5e-12); % Total interval from below
SD.ModFreq = SD.ModFreq / 1e6;

SD.MeasList = genMeasList(SD,'all');

[phi0, Afda] = genBornMat(SD, Medium, [], 'Born', [1 0]);
[phi0, Afds] = genBornMat(SD, Medium, [], 'Born', [0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build a time-domain SD structure and solve the forward problem

disp('Calculating TD problem');

SDtd = SD;
SDtd = rmfield(SDtd,'ModFreq'); % Not relevant

SDtd.TimeDelay     = 5*[0:2047]*1e-12; % Delay time, in sec
SDtd.TimeGateWidth = 0.5*1e-12;        % Time gate width, in sec

SDtd.MeasList = genMeasList(SDtd,'all');

% Call the normal matrix generating routines

% These work...
%[phi1,Atda] = genBornMat(SDtd, Medium, [], 'Born', [1 0]);
%[phi1,Atds] = genBornMat(SDtd, Medium, [], 'Born', [0 1]);

% ...but we can also cut out the middle man here.
[phi1,Atda] = TD3pt(SDtd, Medium, [], [1 0]);
[phi1,Atds] = TD3pt(SDtd, Medium, [], [0 1]);

% Take the interval out of the integrals over the time gate, or the 
% comparison will be off

if (SDtd.TimeGateWidth > 0)
   phi1 = phi1 / SDtd.TimeGateWidth;
   Atda = Atda / SDtd.TimeGateWidth;
   Atds = Atds / SDtd.TimeGateWidth;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fourier transform the time domain and throw up a plot of both amplitude
%  and phase for TD and FD.

clear domus domua

n  = 1:50;
w  = SD.ModFreq;

KT = fft(phi1)*SDtd.TimeDelay(end)/length(phi1);
KF = conj(phi0); % Matlab did something strange?

figure;
subplot(2,3,1); semilogy(w(n)/1000, abs(KF(n)), w(n)/1000, abs(KT(n)))
title('2-pt Green''s function','fontsize',12)
ylabel('Amplitude','fontsize',12);
subplot(2,3,4);   plot(w(n)/1000, angle(KF(n)), w(n)/1000, angle(KT(n)))
axis([w(1)/1000 w(n(end))/1000 -pi pi]);
ylabel('Phase','fontsize',12);
xlabel('Frequency [GHz]','fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KT = fft(Atda)*SDtd.TimeDelay(end)/length(phi1);
KT(1) = conj(KT(1)); % Matlab weirdness
KF = conj(Afda); % Matlab did something strange?

subplot(2,3,2); semilogy(w(n)/1000, abs(KF(n)), w(n)/1000, abs(KT(n)))
title('3-pt Green''s function [abs]','fontsize',12)
ylabel('Amplitude','fontsize',12);
subplot(2,3,5);   plot(w(n)/1000, angle(KF(n)), w(n)/1000, angle(KT(n)))
axis([w(1)/1000 w(n(end))/1000 -pi pi]);
ylabel('Phase','fontsize',12);
xlabel('Frequency [GHz]','fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KT = fft(Atds)*SDtd.TimeDelay(end)/length(phi1);
KT(1) = conj(KT(1)); % Matlab weirdness
KF = conj(Afds); % Matlab did something strange?

subplot(2,3,3); semilogy(w(n)/1000, abs(KF(n)), w(n)/1000, abs(KT(n)))
title('3-pt Green''s function [scat]','fontsize',12)
ylabel('Amplitude','fontsize',12);
subplot(2,3,6);   plot(w(n)/1000, angle(KF(n)), w(n)/1000, angle(KT(n)))
axis([w(1)/1000 w(n(end))/1000 -pi pi]);
ylabel('Phase','fontsize',12);
xlabel('Frequency [GHz]','fontsize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

