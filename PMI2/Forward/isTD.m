% ISTD   Try to guess if this is a Time-domain imager
%
% bool = isTD(SD, MeasList)
%
% SD   PMI structure
% bool = 1 if true, else zero.

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

function[is_td] = isTD(SD, MeasList, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

if (any(any(MeasList(:,6:7) > 0)))
   % Delay time and gate width only make sense in a time-domain context
   
   is_td = 1;

   if (Debug)
     disp('Detected time-domain imager');
   end
   
   is_td = 1;
else
   is_td = 0;
end

return;
