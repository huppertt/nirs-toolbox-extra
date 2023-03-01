% NEWFITBACKGROUND - Iterative search for optical properties
%
% [SD2, Medium2] = newFitBackground(SD, Medium, MeasList, Phi, muvec);
%
% SD, Medium, MeasList - PMI structures
% PhiMeas              - measured data
% muvec                - 2-element vector of flags.  Finds background
%                        absorbtion if muvec(1) is non-zero, background
%                        scattering if muvec(2) is non-zero.
%
% SD2, Medium2         - updated PMI structures with results of fit.

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

function[SD2, Medium2] = newFitBackground(SD, Medium, MeasList, Phi, muvec);

if (isfield(Medium,'Object'))
   error(['Medium.Object cannot be set while finding background optical' ...
	  ' properties' ]);
end

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('muvec','var') | isempty(muvec))
   % Assume you want both Mu's but not N.  There are easier ways to get
   % just one Mu than this.
   
   muvec = [ 1 1 ];
elseif (length(muvec) ~= 2)
   error('muvec must be a 2-element vector');
end

nsrc = size(SD.SrcPos,1);
ndet = size(SD.DetPos,1);
nfrq = length(SD.ModFreq);
nwvl = length(SD.Lambda);

lambda     = 1e-3;
delta_x    = 100;

sdErr      = 0.50;	% Guess the error magnitude
sdErrPhase = 0.50;
muErr      = 0.05;

dyStop     = 1e-5;      % These seem to be reasonable values
lambdaStop = 1e6;
maxLoop    = 250;

Method     = 'Rytov';

% Generate PhiInc for current estimate of SD POS

Phi0 = DPDWHelmholtz(SD, Medium, MeasList);

if (strcmpi(Method, 'Rytov'))
   Y = log(Phi ./ Phi0);
else
   Y = Phi - Phi0;
end

Y2mean     = mean(abs(Y).^2);
Y2mean_old = 1e6*Y2mean;

plotState(SD, Medium, MeasList, Phi, Phi0);

disp(sprintf('Initial y2mean=%g', Y2mean));

% Start with initial values

SD1          = SD;
SD1.MeasList = MeasList;
Medium1      = Medium;
clear SD Medium MeasList;

% Loop until converged or hit maximum

iLoop = 1;

while (abs((Y2mean - Y2mean_old)/Y2mean_old) > dyStop & ...
       lambda < lambdaStop & iLoop < maxLoop)
   
   % Generate Jacobians for this iteration

   muScale = [ Medium1.Muao; Medium1.Muspo ] * muErr;
   
   for k = 1:length(Medium1.Muao)
      tmpml = find(SD1.MeasList(:,4)==k & SD1.MeasList(:,3) > 0);

      if (~isempty(tmpml))
	 fl = unique(SD1.MeasList(tmpml,3));
	 
	 if all(SD1.ModFreq(fl) == 0)
	    % CW measurements are not allowed to change both absorption and
            % scattering, since they're degenerate measures
	    
	    muScale(2,k) = 0;
	 end
      end
      
      clear tmpml fl;
   end
   
   dMu(1,:) = Medium1.Muao / 10;
   dMu(2,:) = Medium1.Muspo/ 10;

   disp('Recalculating Jacobians');
   
   Phi1 = DPDWHelmholtz(SD1, Medium1);

   Jmu0 = genMuJacobian(SD1, Medium1, [], Method, dMu, muvec);
   Jamp = genSDJacobian(SD1, Medium1, [], Phi, Method);

   % Jamp is rank deficient, delete rows to improve the situation.  Delete
   % the first source used for each _unique_ combination of frequency and
   % wavelength.  Leave detectors alone, the solution for unpaired
   % detectors is exp(0.0), which is just what I want.
   
   vFrq = unique(SD1.MeasList(:,3));
   vWvl = unique(SD1.MeasList(:,4));

   nSrc = size(SD1.SrcPos,1);
   nDet = size(SD1.DetPos,1);
   nFrq = length(SD1.ModFreq);
   nWvl = length(SD1.Lambda);
   
   ns = length(unique(SD1.MeasList(:,1)));
   nd = length(unique(SD1.MeasList(:,2)));

   if (ns > 1 & nd > 1)
      % Multiple sources and multiple detectors
      
      for iF = 1:length(vFrq)
	 iFrq = vFrq(iF);
	 
	 for iW = 1:length(vWvl)
	    iWvl = vWvl(iW);
	    
	    tmpml = find(SD1.MeasList(:,3) == iFrq & ...
			 SD1.MeasList(:,4) == iWvl);
	    iSrc = SD1.MeasList(tmpml(1),1);
	    
	    iJS = (iFrq-1)*nSrc*nWvl + (iWvl-1)*nSrc + iSrc;
	 
	    Jamp(tmpml, iJS) = 0;
	 end
      end
   else
      % Exactly ones source and/or one detector - no hope
      error('Cannot do fitting with only one source and/or one detector');
   end
   
   clear iFrq iWvl iSrc iJS iJD tmpml;
   
   % Invert Jacobians to provide new best guess

   [SD2, Medium2, delta_x] = ...
       UpdateX(Jmu0, Jamp, Y, lambda, SD1, Medium1, [], ...
	       sdErr, sdErrPhase, muScale, muvec);
   
   % Calculate residue (goodness-of-fit)
   
   if (all(Medium2.Muao > 0) & all(Medium2.Muspo > 0))
      Phi1 = DPDWHelmholtz(SD2, Medium2);
      
      if (strcmpi(Method, 'Rytov'))
	 Ynew = log(Phi ./ Phi1);
      else
	 Ynew = Phi - Phi1;
      end
      
      Y2mean_old = Y2mean;
      Y2mean     = mean(abs(Ynew).^2);
   end   

   % Accept/Reject the move
   
   if ((Y2mean < Y2mean_old) & all(Medium2.Muao > 0) & all(Medium2.Muspo > 0))
      % Accept the move
      
      fprintf(1, 'y2mean=%f lambda=%g Musp=[ ', Y2mean, lambda);
      fprintf(1, '%5.2f ', Medium2.Muspo);
      fprintf(1, '] Mua=[ ');
      fprintf(1, '%5.3f ', Medium2.Muao);
      fprintf(1, ']\n');

      lambda = lambda / 2;
      Y      = Ynew;
      iLoop  = iLoop + 1;
   else
      while ((lambda < lambdaStop) & ...
	     ((Y2mean > Y2mean_old) | ...
	      any(Medium2.Muao <= 0) | any(Medium2.Muspo <= 0)))

	 % Reject the move, try again
	 
	 fprintf(1, '---y2mean=%f lambda=%g Musp=[ ', Y2mean, lambda);
	 fprintf(1, '%5.2f ', Medium2.Muspo);
	 fprintf(1, '] Mua=[ ');
	 fprintf(1, '%5.3f', Medium2.Muao);
	 fprintf(1, ']\n');

	 lambda = lambda * 10;
	 
	 [SD2, Medium2, delta_x] = ...
             UpdateX(Jmu0, Jamp, Y, lambda, SD1, Medium1, [], ...
		     sdErr, sdErrPhase, muScale, muvec);

	 if (all(Medium2.Muao > 0) & all(Medium2.Muspo > 0))
	    Phi1 = DPDWHelmholtz(SD2, Medium2);
	    
	    if (strcmpi(Method, 'Rytov'))
	       Ynew = log(Phi ./ Phi1);
	    else 
	       Ynew = Phi - Phi1;
	    end
	    
	    Y2mean = mean(abs(Ynew).^2);

	    % Accept the move
	    
	    if ((Y2mean < Y2mean_old) & ...
		all(Medium2.Muao > 0) & all(Medium2.Muspo > 0))
	       
	       fprintf(1, '+++y2mean=%f lambda=%g Musp=[ ', Y2mean, lambda);
	       fprintf(1, '%5.2f ', Medium2.Muspo);
	       fprintf(1, '] Mua=[ ');
	       fprintf(1, '%5.3f', Medium2.Muao);
	       fprintf(1, ']\n');

	       lambda = lambda / 2;
	       Y      = Ynew;
	       iLoop  = iLoop + 1;
	    end
	 else
	    % Rejected step, preserve previous ymean
	    Y2mean = Y2mean_old;
	 end
      end   
      
      if (lambda >= lambdaStop & iLoop > 1)
	 % Bad estimate, go with previous guess
	 
	 SD2     = SD1;
	 Medium2 = Medium1;
	 break;
      end
   end

   SD1     = SD2;
   Medium1 = Medium2;
   clear SD2 Medium2;
   
   disp(sprintf('delta_x -> %g, dy = %g', ...
		max(abs(delta_x(:))), abs((Y2mean-Y2mean_old)/ Y2mean_old)));

   plotState(SD1, Medium1, [], Phi, Phi1);
end

% Pass final results back to caller

SD2     = SD1;
Medium2 = Medium1;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do one update cycle using the parameter "lambda" for regularization
%

function[SD2, Medium2, delta_x] = ...
    UpdateX(Jmu0, Jamp, Y, lambda, SD1, Medium1, MeasList, ...
            ampErr, ampErrPhase, muScale, muvec)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD1.MeasList;
end

nSrc = size(SD1.SrcPos,1);
nDet = size(SD1.DetPos,1);
nWvl = max(1,length(SD1.Lambda));
nFrq = max(1,length(SD1.ModFreq));

% Pre-condition matrices

Jmu0 = Jmu0 * diag([ muScale(1,:) muScale(2,:) ]);
Jamp = Jamp * ampErr;

% Pack individual Jacobians into a single large array

if (all(SD1.ModFreq == 0))
   J = [ Jmu0 Jamp ];
else
   J = [ Jmu0 Jamp i*Jamp*(ampErrPhase/ampErr) ];
   
   J = [ real(J); imag(J) ];
   Y = [ real(Y); imag(Y) ];
end

% Calculate the update using non-singular parts of the matrix

disp('Inverting combined Jacobian');

nzel = find(sum(J.^2) > 0);

if (1)
   if (size(J,1) > length(nzel))
      % More measurements than unknowns

      delta_x0 = tik(J(:,nzel), Y, lambda);
   else
      % More unknowns than measurements
      
      delta_x0 = tik(J(:,nzel), Y, lambda);
   end
   
   %  B = B + lambda * max(diag(B)) * speye(size(B));
   % delta_x0 = B \ (J(:,nzel)' * Y);
else
   % This should work just as well and be MUCH faster
   
   Jnorm = normest(J); % Magnitude of largest eigenvalue of J
   
   B = [ J(:,nzel); lambda * Jnorm * speye(length(nzel)) ];
   
   Yp = Y; 
   Yp(end+1:size(B,1)) = 0;

   delta_x0 = tcgls(B, Yp, [], 100);
end

disp('Repacking update');

% Copy back into original ordering

delta_x = zeros(1,size(J,2));
delta_x(nzel) = delta_x0;

% Initialze new SD/Medium with old values

SD2     = SD1;
Medium2 = Medium1;
nCol    = 0;

SD2.MeasList = MeasList;

dmu0 = delta_x(nCol + [1:size(Jmu0,2)]) * diag([ muScale(1,:) muScale(2,:) ]);

Medium2.Muao  = Medium1.Muao  + dmu0(1:end/2);
Medium2.Muspo = Medium1.Muspo + dmu0(end/2+1:end);
nCol = nCol + size(Jmu0,2);

SD2.SrcAmp = SD1.SrcAmp .* ...
    exp(reshape(delta_x(nCol + [1:nSrc*nWvl*nFrq]), nSrc, nWvl, nFrq)*ampErr);
nCol = nCol + nSrc*nWvl*nFrq;

SD2.DetAmp = SD1.DetAmp .* ...
    exp(reshape(delta_x(nCol + [1:nDet*nWvl*nFrq]), nDet, nWvl, nFrq)*ampErr);
nCol = nCol + nDet*nWvl*nFrq;

if (any(SD1.ModFreq) ~= 0)
   SD2.SrcAmp = SD2.SrcAmp .* ...
       exp(i * reshape(delta_x(nCol + [1:nSrc*nWvl*nFrq]), ...
		       nSrc, nWvl, nFrq) * ampErrPhase);
   nCol = nCol + nSrc*nWvl*nFrq;

   SD2.DetAmp = SD2.DetAmp .* ...
       exp(i * reshape(delta_x(nCol + [1:nDet*nWvl*nFrq]), ...
		       nDet, nWvl, nFrq) * ampErrPhase);
   nCol = nCol + nDet*nWvl*nFrq;
end

% We can still turn purely real-valued SD's into complex numbers due
% to round-off errors in the inversion.  Fix this here.

fl = find(SD1.ModFreq == 0);

if (~isempty(fl))
   SD2.SrcAmp(:,:,fl) = real(SD2.SrcAmp(:,:,fl));
   SD2.DetAmp(:,:,fl) = real(SD2.DetAmp(:,:,fl));
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% internal function to plot the progress
%
function plotState(SD, Medium, MeasList, PhiExp, PhiInc)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

figure(1);

% plot the fluence versus separation

dR = calcSep(SD, MeasList);
m1 = floor(min(dR));
m2 = ceil(max(dR));

sdCor = zeros(size(MeasList,1),1);

for k = 1:size(MeasList,1)
   si = MeasList(k,1);
   di = MeasList(k,2);
   fi = MeasList(k,3);
   wi = MeasList(k,4);
   
   sdCor(k) = SD.SrcAmp(si,wi,fi) .* SD.DetAmp(di,wi,fi);
end

if (SD.ModFreq == 0)
   subplot(1,1,1)
   plot(dR, log(PhiExp./PhiInc),'r.')
   set(gca,'xlim',[m1 m2]);
else
   subplot(1,2,1)
   plot(dR, log(abs(PhiExp./PhiInc)), 'r.');

   lh = line([m1 m2],[0 0]);
   set(lh,'color','black','linestyle','--');
   set(gca,'xlim',[m1 m2],'ylim',[-1 1]*max(abs(get(gca,'ylim'))));
   
   subplot(1,2,2)
   plot(dR, angle(PhiExp./PhiInc), 'r.');

   lh = line([m1 m2],[0 0]);
   set(lh,'color','black','linestyle','--');
   set(gca,'xlim',[m1 m2],'ylim',[-1 1]*max(abs(get(gca,'ylim'))))
end

drawnow

return
