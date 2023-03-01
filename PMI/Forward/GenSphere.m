%GENSPHERE      Create an absorbing sphere in a discrete volume.
%
%   x = gensphere(CompVol, ctr, radius, value)
%
%   x           The object function.
%
%   CompVol     A structure defining the computational volume.  This
%               structure should have the members: Type, X, Y and Z.  Type
%               should be uniform specifying a uniform sampling volume of
%               voxels. X, Y and Z are vectors specifying the centers of the
%               voxels.
%
%   ctr         The center of the sphere.
%
%   radius      The radius of the sphere.
%
%   value       The value within the boundaries of the sphere.

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: GenSphere.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 17:39:56  rjg
%  Initial Revision for PMI 3.0
%
%  Revision 2.0  1998/08/05 15:50:55  rjg
%  Start of version 2
%  Handles structure for computational volume.
%
%  Revision 1.1  1998/06/03 16:03:58  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = GenSphere(CompVol, ctr, radius, value)

%%
%%  Create a disk centered set of coordinates
%%
if strcmp(CompVol.Type, 'uniform') == 0
    error('This routine only handles uniform voxelations');
end
[X Y Z] = meshgrid(CompVol.X-ctr(1), CompVol.Y-ctr(2), CompVol.Z-ctr(3));

%%
%%  Compute the spherical radius of each coordinate and find all that are ...
%%  within the radius and set them to the delta mu_a
%%
r = sqrt(X.^2 + Y.^2 + Z.^2);
X =  (r <= radius) * value;

