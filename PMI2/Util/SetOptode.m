% SETOPTODE   Set a grid of optode positions and amplitudes.
%             Optodes will be placed at all combinations of X, Y,
%             and Z.
%
%   [Pos, Amp] = SetOptode(X, Y, Amp)
%
%   X, Y - vectors of positions [in cm] for the appropriate axis
%
%   Amp - is a matrix of dimension [#Optodes, #Wavelengths] that 
%         specifies the wavelength dependent amplitude of each
%         optode.  If this is of dimension [1, #Wavelengths], then
%         the same amplitude will be given to all optodes for each
%         wavelength.  The units for the source optode is Watts
%         entering the tissue.  The units for the detector is the
%         detector aperture in cm^2 times the coupling efficiency
%         of the detector.
%
% Returns:
%   Pos - Matrix of dimension [#Optodes, 3] for use as SD.SrcPos
%         or SD.DetPos
%
%   Amp - is a matrix of dimention [#Optodes, #Wavelengths, 1] giving
%         the wavelength dependent amplitude of each optode.  MAY NOT
%         WORK CORRECTLY IF YOU HAVE MULTIPLE FREQUENCIES!

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

function [Pos, Amplitude] = SetOptode(X, Y, Z, Amp);

[Xp Yp Zp] = meshgrid(X, Y, Z);
Pos = [Xp(:) Yp(:) Zp(:)];
num = size(Pos,1);
  
if size(Amp,1) == num
   Amplitude = Amp;
else
   Amplitude = ones(num,1) * Amp;
end

return;
