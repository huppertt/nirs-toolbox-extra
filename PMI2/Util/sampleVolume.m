% SAMPLEVOLUME  Compute sample volume from Medium.CompVol
%
% [X, Y, Z, dV] = sampleVolume(Medium.CompVol);
%
% Medium - PMI Structure
%
% X,Y,Z  - coordinates of center of each voxel
% dV     - volume of the voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Xm, Ym, Zm, dV] = sampleVolume(CompVol)

if (~isstruct(CompVol))
   error('Medium.CompVol is a structure');
end

% Based on CompVol.Type, generate the sample volume and the (uniform!)
% volume of each voxel

if     strcmpi(CompVol.Type, 'uniform')
   [Xm Ym Zm] = meshgrid(CompVol.X, CompVol.Y, CompVol.Z);
   
   dV = CompVol.XStep * CompVol.YStep * CompVol.ZStep;
elseif strcmpi(CompVol.Type, 'computed')
   dX = CompVol.XStep;
   dY = CompVol.YStep;
   dZ = CompVol.ZStep;
   
   % The range can be specified either way

   if (isfield(CompVol,'X'))
      X = [CompVol.X(1) : dX : CompVol.X(2)];
   else
      X = [CompVol.X1 : dX : CompVol.X2];
   end

   if (isfield(CompVol,'X'))
      Y = [CompVol.Y(1) : dY : CompVol.Y(2)];
   else
      Y = [CompVol.Y1 : dY : CompVol.Y2];
   end

   if (isfield(CompVol,'X'))
      Z = [CompVol.Z(1) : dZ : CompVol.Z(2)];
   else
      Z = [CompVol.Z1 : dZ : CompVol.Z2];
   end
   
   [Xm Ym Zm] = meshgrid(X, Y, Z);
   dV = dX * dY * dZ;
   
   clear dX dY dZ X Y Z;
elseif strcmpi(CompVol.Type, 'list')
   Xm = CompVol.X;
   Ym = CompVol.Y;
   Zm = CompVol.Z;
   
   dV = CompVol.Volume;
else
   error(['Unknown volume specification ' CompVol.Type ]);
end

Xm = Xm(:);
Ym = Ym(:);
Zm = Zm(:);

return;

