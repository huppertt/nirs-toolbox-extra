% CALCSEP  Calculate separation between optode pairs in measurement list.  
%
% A very simple calculation, but one that is used often.
%
% dR = calcSep(SD, [MeasList]);

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

function[dR] = calcSep(SD, MeasList)

if (~exist('MeasList','var') | isempty(MeasList))
  % use default
  MeasList = SD.MeasList;
end

sl = MeasList(:,1);
dl = MeasList(:,2);

dR = sqrt(sum((SD.SrcPos(sl,:) - SD.DetPos(dl,:)).^2, 2));

return;
