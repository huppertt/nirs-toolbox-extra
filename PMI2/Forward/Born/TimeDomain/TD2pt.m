% TD2PT  Generate the incident field for TD model
%
% Phi0 = TD2pt(SD, Medium, Measlist, [Debug]);
%   
%   SD, Medium - PMI data structures
%   MeasList   - Measurements to use
%
% Returns:
%    Phi0 - The total measured fluence for each
%         source-detector pair specified in SD.MeasList given a
%         homogneeous medium.  

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

function[Phi0] = TD2pt(SD, Medium, MeasList, Debug)

if (~exist('Debug','var'))
   Debug = 0;
end

if ((~exist('MeasList','var')) | (isempty(MeasList)))
   MeasList = SD.MeasList;
end

if (length(Medium.Muao) < length(unique(MeasList(:,4))))
   % This is a common error, check for it explicitly
   error('More wavelengths than optical properties');
end

switch lower(Medium.Geometry)
   case {'infinite', 'inf'}
      if (exist('mextd2ptNB','file') == 3)
         Phi0 = mextd2ptNB(SD, Medium, MeasList, Debug);
      else
         Phi0 = TD_2ptNB(SD, Medium, MeasList, Debug);
      end

   case { 'semi-infinite', 'semi', 'extrapolated'}
      if (exist('mextd2ptZB','file') == 3)
         Phi0 = mextd2ptZB(SD, Medium, MeasList, Debug);
      else
         Phi0 = TD_2ptZB(SD, Medium, MeasList, Debug);
      end

   case {'slab'}
      if (exist('mextd2ptSB','file') == 3)
         Phi0 = mextd2ptSB(SD, Medium, MeasList, Debug);
      else
         Phi0 = TD_2ptSB(SD, Medium, MeasList, Debug);
      end

   otherwise
      error(['Unknown or unsupported boundary condition: ' Medium.Geometry]);
end

return;
