% FITBACKGROUND - brute force search for optical properties
%
% [mus, mua, cf] = fitBackground(SD, Medium, MeasList, PhiMeas);
%
% SD, Medium, MeasList - PMI structures
% PhiMeas              - measured data

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

function[mus, mua, cost] = fitBackground(SD, Medium, MeasList, PhiMeas, ...
                                         doAmp, doPlot, Debug);

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('doAmp','var') | isempty(doAmp))
   % Do amplitude fitting 
   doAmp = 1;
end

if (~exist('doPlot','var') | isempty(doPlot))
   % Status plot
   doPlot = 1;
end

if (~exist('Debug','var') | isempty(Debug))
   % Extra debugging info
   Debug = 0;
end

% Search ranges

mua = [ 0.01:0.01:0.15 ];
mus = [ 2:20 ];

if (isfield(SD,'ModFreq') & (SD.ModFreq ~= 0))
   hasPhs = 1;
else
   hasPhs = 0;
end

% Initialize tables

if (hasPhs)
   cost = zeros(length(mua), length(mus), length(SD.Lambda), 2);
else
   cost = zeros(length(mua), length(mus), length(SD.Lambda));
end
      
% Do each wavelength separately

%for l = 1:length(SD.Lambda)
%   figure(l);
%   axis([ mus(1) mus(end) mua(1) mua(end) ]);
%   set(l,'DoubleBuffer','on');
%end

if (doPlot)
   ploth = figure;
   set(ploth,'doublebuffer','on');
else
   ploth = [];
end

% List of figures used
figh = zeros(1,length(SD.Lambda));

for ai = 1:length(mua)
   Medium.Muao  = ones(1,length(SD.Lambda)) * mua(ai);

   for si = 1:length(mus)
      Medium.Muspo = ones(1,length(SD.Lambda)) * mus(si);

      if (Debug)
         disp(sprintf('Checking [mua,mus] = (%f,%f)', mua(ai), mus(si)));
      end

      k = (ai-1) * length(mus) + si;

      for l = 1:length(SD.Lambda)
	 ml = find(MeasList(:,4)==l);
	 
         if (~isempty(ml))
	    tmpML = MeasList(ml,:);
 
	    c = CheckMu(SD, Medium, tmpML, PhiMeas(ml), ...
                        doAmp, doPlot, ploth, Debug);

            for m = 1:length(c)
	       cost(ai,si,l,m) = c(m);
            end

            if (doPlot)
               if (figh(l) <= 0)
                  figh(l) = figure;
		  set(figh(l), 'doublebuffer','on');
               else
                  figure(figh(l));
               end

               imagesc(mus, mua, cost(:,:,l,1));
               drawnow;
            end
         end
      end
   end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check a single [mus,mua] pair, return Chisqr and some debugging terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cost] = CheckMu(SD, Medium, MeasList, Phi, ...
                          doamp, doplot, ploth, debug)

if (isempty(Phi))
  warning('Empty data vector');
  keyboad;
end

if (isfield(SD,'ModFreq') & (any(SD.ModFreq(unique(MeasList(:,3))) ~= 0)))
   hasPhs = 1;
else
   hasPhs = 0;
end

% Get incident fluence, don't need full forward matrix for chisqr

Phi0 = FD2pt(SD, Medium, MeasList);

if ((hasPhs == 0) & min(Phi0) < 0)
   % Negative real/imaginary allowed with RF data, otherwise Phi0>=0
   warning('genBornData() returned a negative fluence');
   keyboard;
end

% Rescale the (possibly corrected) incident fluence

if (doamp)
   [sA, dA, pA] = fitSD(SD, Medium, MeasList, Phi, Phi0);
   Phi0         = Phi0 .* pA;
   
   clear sA dA pA;
else
   Phi0         = Phi0 .* abs(mean(Phi)) ./ abs(mean(Phi0));
   
   if (~isreal(Phi0))
      Phi0 = Phi0 .* exp(i*mean(angle(Phi0./Phi)));
   end
end

% Compute cost function for this set of parameters

cost(1) = mean(log(abs(Phi./Phi0)).^2);

if (hasPhs)
   cost(2) = mean(angle(Phi./Phi0).^2);
end

if (debug)
   disp(sprintf('  ..chisqr = %e', cost));
end

if (doplot)
   h = gcf;
   figure(ploth);
   dR = calcSep(SD,MeasList);
   n = length(Phi);
   subplot(1,2,1); semilogy(1:n, abs(Phi), 'o', 1:n, abs(Phi0));
   subplot(1,2,2); plot(1:n, angle(Phi), 'o', 1:n, angle(Phi0));
                   axis([1 n -pi pi]);
   drawnow;
   figure(h);
end

return;
