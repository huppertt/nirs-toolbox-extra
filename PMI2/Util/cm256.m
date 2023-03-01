% CM256  256 color version of the default "jet" colormap
%
% colormap(cm256());

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

function[cm] = cm256()

cm = zeros(256,3);

cm(  1:97, 3) = min(1, ([1:97]' + 35)/64);
cm( 98:256,3) = max(0, (64 - [1:159]')/64);

cm(  1:92, 2) = max(0, (64 - [92:-1:1]')/64);
cm( 93:256,2) = max(0, min(1,  (225-[93:256]')/64));

cm(  1:157,1) = max(0, ([1:157]'  -  93)/64);
cm(158:256,1) = min(1, (161+128 - [158:256]')/64);

return;