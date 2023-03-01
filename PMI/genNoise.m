%GenNoise       Generate the noise vector and noise scattered field.
%
%   pmi = gennoise(pmi);
%
%   pmi         The Photon Migration Imaging data structure to updated.
%
%
%   GenNoise generates instances the selected noise models in the PMI data
%   structure.
%
%   Calls: none.
%
%   Bugs: none known.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%  $Author: dboas $
%
%  $Date: 2000/09/01 19:32:58 $
%
%  $Revision: 1.5 $
%
%  $Log: genNoise.m,v $
%  Revision 1.5  2000/09/01 19:32:58  dboas
%  Having some trouble with re-initializing the P sub-structure.  This is probably
%  still not the desired fix.
%
%  Revision 1.4  2000/08/16 17:57:27  dboas
%  We don't want to copy ds.Fwd.P.PhiInc to ds.Inv.P.PhiInc.  We only copy
%  PhiTotal.  Also, PhiTotalN always needs to be initialized to PhiTotal before
%  adding noise.
%
%  Revision 1.3  2000/08/02 20:17:36  dboas
%  Initialize the TotalVar variable to ones if there is no noise.
%
%  Revision 1.2  2000/07/27 14:57:20  dboas
%  Add P(idxLambda) for the MODEL and NOISE.
%  Removed a number of noise types that aren't relevant.  We are only
%     keeping ShotNoise and ElectronicNoise.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.5  2000/04/18 19:48:29  dboas
%
%  Small changes to ElectronicSNR and ShotSNR
%
%  Revision 3.4  2000/04/14 22:08:47  rjg
%  Cleared any previous results stored in Noise vectors and STD.
%
%  Revision 3.3  2000/04/13 21:06:10  dboas
%
%  Added two new noise types: Shot noise and Electronic noise.  These are
%  much more realistic noise models for our simulations.  Shot Noise is activated
%  by .ShotSNRflag and controlled by .ShotB which is the detection bandwidth
%  in Hz.  Electronic Noise is activated by .ElectronicSNRflag and controlled by
%  .ElectronicNoise which is the noise equivalent power in Watts.
%
%  Revision 3.2  2000/03/03 04:26:21  rjg
%  Added calculation of .Inv.PhiScat
%
%  Revision 3.1  1999/12/22 16:00:57  rjg
%  Added computation of .Inv.PhiScatw
%
%  Revision 3.0  1999/06/17 19:29:38  rjg
%  Initial Revision for PMI 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmi = genNoise(pmi);

nLambda = length(pmi.Fwd.Lambda);

%%
%%  Clear any previous noise results
%%
idx = 1;

if isfield(pmi.Noise, 'P')
  if isfield(pmi.Noise.P(idx), 'ShotNoiseSTD')
    pmi.Noise.P(idx) = rmfield(pmi.Noise.P(idx), 'ShotNoiseSTD');
  end
  if isfield(pmi.Noise.P(idx), 'ShotNoise')
    pmi.Noise.P(idx) = rmfield(pmi.Noise.P(idx), 'ShotNoise');
  end
  if isfield(pmi.Noise.P(idx), 'ElectronicNoiseSTD')
    pmi.Noise.P(idx) = rmfield(pmi.Noise.P(idx), 'ElectronicNoiseSTD');
  end
  if isfield(pmi.Noise.P(idx), 'ElectronicNoise')
    pmi.Noise.P(idx) = rmfield(pmi.Noise.P(idx), 'ElectronicNoise');
  end
  if isfield(pmi.Noise.P(idx), 'TotalVar')
%    foo = rmfield(pmi.Noise.P(idx), 'TotalVar');
    pmi.Noise.P(idx).TotalVar = [];
  end
  if isfield(pmi.Noise.P(idx), 'w')
%    pmi.Noise.P(idx) = rmfield(pmi.Noise.P(idx), 'w');
    pmi.Noise.P(idx).w = [];
  end
end

%%
%%  Copy in the non-noise scattered field into the noise modified
%%  scattered field.
%%
for idx=1:nLambda
%  pmi.Inv.P(idx).PhiInc(:,1) = pmi.Fwd.P(idx).PhiInc(:,1); THIS IS 
%  DONE BY genFwdMat(ds.Inv)
  pmi.Inv.P(idx).PhiTotal(:,1) = pmi.Fwd.P(idx).PhiTotal(:,1);
  P(idx).TotalVar = zeros(size(pmi.Fwd.P(idx).PhiTotal));
  P(idx).nMeas = size(pmi.Fwd.P(idx).PhiTotal);
end

%if ~pmi.Noise.ShotSNRflag & ~pmi.Noise.ElectronicSNRflag
% THIS LOOP SHOULD ALWAYS BE DONE... RIGHT?
for idx=1:nLambda
  pmi.Inv.P(idx).PhiTotalN(:,1) = pmi.Inv.P(idx).PhiTotal(:,1);
end;
%end;


%%
%%  Add shot noise
%%
if pmi.Noise.ShotSNRflag
    planckc = 6.6e-34 * 3e17;
    for i=1:nLambda
      ShotNoiseSTD  = planckc/pmi.Fwd.Lambda(i) .* ...
          (pmi.Noise.ShotSNR * abs(pmi.Fwd.P(i).PhiTotal) .* ...
           pmi.Fwd.Lambda(i)/planckc).^0.5;

      P(i).TotalVar = P(i).TotalVar + ShotNoiseSTD .^ 2;
      ShotNoise = randn(P(i).nMeas) .* ShotNoiseSTD;
      pmi.Inv.P(i).PhiTotalN = pmi.Inv.P(i).PhiTotalN + ShotNoise;
      if pmi.Debug
        pmi.Noise.P(i).ShotNoiseSTD = ShotNoiseSTD;
        pmi.Noise.P(i).ShotNoise = ShotNoise;
        fprintf('Max shot noise std %e\n', max(pmi.Noise.P(i).ShotNoiseSTD));
        fprintf('Min shot noise std %e\n\n', min(pmi.Noise.P(i).ShotNoiseSTD));
      end
      
      if length(find(pmi.Inv.P(i).PhiTotal./ShotNoiseSTD)<1)>0
	disp('WARNING: The SNR due to shot noise is less than 1 for some measurements');
      end
    end
end

%%
%%  Add detector electronic noise
%%
if pmi.Noise.ElectronicSNRflag
  for idx=1:nLambda
    ElectronicNoiseSTD = ones(size(pmi.Fwd.P(idx).PhiTotal)) * pmi.Noise.ElectronicSNR;
    P(idx).TotalVar = P(idx).TotalVar + ElectronicNoiseSTD .^ 2;
    ElectronicNoise = randn(P(idx).nMeas) .* ElectronicNoiseSTD;
    pmi.Inv.P(idx).PhiTotalN = pmi.Inv.P(idx).PhiTotalN + ElectronicNoise;
    if pmi.Debug
      pmi.Noise.P(idx).ElectronicNoiseSTD = ElectronicNoiseSTD;
      pmi.Noise.P(idx).ElectronicNoise = ElectronicNoise;
      fprintf('Max electronic noise std %e\n', max(pmi.Noise.P(idx).ElectronicNoiseSTD));
      fprintf('Min electronic noise std %e\n\n', min(pmi.Noise.P(idx).ElectronicNoiseSTD));
    end

    if length(find(pmi.Inv.P(idx).PhiTotal./ElectronicNoiseSTD)<1)>0
      disp('WARNING: The SNR due to electronic noise is less than 1 for some measurements');
    end
  end
end


%
% This information is now needed in procNoise.  DAB 2000-05-01
%
%if pmi.Debug
for idx = 1:nLambda
    if P(idx).TotalVar(1) == 0
      pmi.Noise.P(idx).TotalVar= ...
	  ones(size(pmi.Fwd.P(idx).PhiTotal));
    else
      pmi.Noise.P(idx).TotalVar = P(idx).TotalVar;
    end
end
%end


