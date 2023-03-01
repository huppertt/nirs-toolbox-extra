% GENBLOCK   Create an absorbing block in a discrete volume.
%
% x = genBlock(Medium, ctr, dim, value)
%
%   Medium      PMI structure ...
%
%   ctr         The center of the block as [cX cY cZ].
%
%   dim         The length of each side of the block as [lX lY lZ].
%
%   value       The value within the boundaries of the sphere.
%
% Returns:
%   x           The object function.

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

function[X] = genBlock(Medium, ctr, dim, value)

%%
%%  Create a block centered set of coordinates
%%

[X Y Z dV] = sampleVolume(Medium.CompVol);

X = X - ctr(1);
Y = Y - ctr(2);
Z = Z - ctr(3);

%%
%%  Compute the spherical radius of each coordinate and find all that are ...
%%  within the radius and set them to the delta mu_a
%%

BlkX = logical(abs(X) <= dim(1)/2);
BlkY = logical(abs(Y) <= dim(2)/2);
BlkZ = logical(abs(Z) <= dim(3)/2);

X = logical(BlkX & BlkY & BlkZ) * value;
X = X(:);

return;
