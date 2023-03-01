% FD2PT           Generate the incident field for FD model
%
% Phi0 = FD2pt(SD, Medium, Measlist, [Debug]);
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

function[Phi0] = FD2pt(SD, Medium, MeasList, Debug)

if ((~exist('MeasList','var')) | (isempty(MeasList)))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var'))
   Debug = 0;
end

% Call the appropriate forward method based on the geometry

switch lower(Medium.Geometry)
   case { 'infinite', 'inf', 'inft' }
      if Debug
	 fprintf(['Executing infinite medium boundary computation\n']);
      end
      
      Phi0 = Hlm2ptNB(SD, Medium, MeasList, Debug);
      
   case { 'semi-infinite', 'semi', 'extrapolated' }
      if Debug
	 fprintf(['Executing extrapolated zero boundary computation\n']);
      end

      Phi0 = Hlm2ptZB(SD, Medium, MeasList, Debug);
      
   case { 'slab' }
      if Debug
	 fprintf(['Executing slab boundary computation\n']);
      end

      Phi0 = Hlm2ptSB(SD, Medium, MeasList, Debug);
      
   otherwise
      error(['Unknown boundary condition: ' Medium.Geometry]);
end

return;


