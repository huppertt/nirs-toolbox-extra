% MEXTD3PTZB 1st Born approx. of a perturbed Helmholtz equation.
%
%   MEXTD3PTZB computes the forward weighting matrix associated
%   with the Born-1 approximation to a spatial varying k and semi-
%   infinite boundary conditions.
%
%   Do not specify the sources or detector exactly at a sampling point to
%   prevent dived by zero errors when evaluating the Green's functions.
%
% A = mextd3ptZB(SD, Medium, MeasList, OptProp, Debug);
%
% SD, Medium - PMI data structures
% MeasList   - Measurement list or [] to use SD.MeasList
% OptProp    - Optical properties to model
%              OptProp(1) ~= 0 -> Absorping perturbations
%              OptProp(2) ~= 0 -> Scattering perturbations
%
% Returns:
%   A        The forward matrix relating the contribution from each
%            voxel to a specific source-detector pair.  Each column
%            is for a different measurement.  The perturbation
%            is \delta{\mu} (abs) or \frac{\delta D}{D} (scat). 

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
