% MLtoFB  Given a measurement list, calculate the full-born indicies.
%         Common code, to be shared by the different full-born routines.
%
% FB = MLtoFB(SD, Medium, MeasList, OptProp);
%
% FB - structure with parameters for mapping MeasList to FullBorn indices
%
% DO NOT CALL THIS CODE DIRECTLY!

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

function [FB] = MLtoFB(SD,Medium,MeasList,OptProp)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

% Scattering vs absorption

FB.calc_mua  = OptProp(1);
FB.calc_musp = OptProp(2);

%  Create the sampling volume

[Xm, Ym, Zm, volVoxel] = sampleVolume(Medium.CompVol);

nVox = length(Xm);
rVox = [ Xm Ym Zm ];   % Nx3 vector of voxel positions
R    = (volVoxel / (4*pi/3))^(1/3);

FB.nVox = nVox;
FB.rVox = rVox;
FB.voxR = R;
FB.volVoxel = volVoxel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find all the unique combinations of optical properties in the
%  measurement list (independant of source/detector index)

MLtmp = MeasList; 
MLtmp(:,[1 2]) = 0;

% Number of unique sets of optical properties (Green's functions)
[MLopt,optF,optR] = unique(MLtmp,'rows');
nOpt  = size(MLopt, 1);
nPts  = nOpt * nVox;

FB.MLopt = MLopt;
FB.optF  = optF;
FB.optR  = optR;
FB.nOpt  = nOpt;
FB.nPts  = nOpt * nVox;

clear MLtmp

% Find all the unique source combinations

MLtmp = MeasList; MLtmp(:,2) = 0;

[MLsrc,srcF,srcR] = unique(MLtmp,'rows');
nSrcs = size(MLsrc,1);

FB.MLsrc = MLsrc;
FB.srcF  = srcF;
FB.srcR  = srcR;
FB.nSrcs = nSrcs;

clear MLtmp

% Find all the unique detector combinations

MLtmp = MeasList; MLtmp(:,1) = 0;

[MLdet,detF,detR] = unique(MLtmp,'rows');
nDets = size(MLdet,1);

FB.MLdet = MLdet;
FB.detF  = detF;
FB.detR  = detR;
FB.nDets = nDets;

clear MLtmp

return;
