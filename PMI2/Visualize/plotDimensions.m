% PLOTDIMENSIONS  Add some aesthetics to figures with many subplots.  
%
% [M,N] = plotDimensions(nz, doMua, doMus);
%
% Not normally called from user code.  For large nz (nz >= 25), falls
%  M,N are roughly the sqrt(nz), rounded up if needed.

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

function[M,N] = plotDimensions(nz, doMua, doMus)

if (doMua & doMus)
   % If BOTH doMua and doMus are defined then M must be even
   % Cutoff here is nz=18 for each instead of nz=25
   
   switch (nz)
      case 1
	 M = 2; N = 1;
      case 2
	 M = 2; N = 2;
      case 3
	 M = 2; N = 3;
      case 4
	 M = 2; N = 4;
      case 5
	 M = 4; N = 3;
      case 6
	 M = 4; N = 3;
      case 7
	 M = 4; N = 4;
      case 8
	 M = 4; N = 4;
      case 9
	 M = 4; N = 5;
      case 10
	 M = 4; N = 5;
      case 11
	 M = 4; N = 6;
      case 12
	 M = 4; N = 6;
      case 13
	 M = 4; N = 7;
      case 14
	 M = 4; N = 7;
      case 15
	 M = 6; N = 5;
      case 16
	 M = 6; N = 6;
      case 17
	 M = 6; N = 6;
      case 18
	 M = 6; N = 6;
      otherwise
	 M = 2*ceil(sqrt(nz)); N = ceil(sqrt(nz));
   end
else
   % If EITHER doMua or doMus are defend, M may be even or odd

   switch (nz)
      case 1
	 M = 1; N = 1;
      case 2
	 M = 1; N = 2;
      case 3
	 M = 2; N = 2;
      case 4
	 M = 2; N = 2;
      case 5
	 M = 2; N = 3;
      case 6
	 M = 2; N = 3;
      case 7
	 M = 3; N = 3;
      case 8
	 M = 3; N = 3;
      case 9
	 M = 3; N = 3;
      case 10
	 M = 3; N = 4;
      case 11
	 M = 3; N = 4;
      case 12
	 M = 3; N = 4;
      case { 13, 14, 15, 16 }
	 M = 4; N = 4;
      case { 17, 18, 19, 20 }
	 M = 4; N = 5;
      case { 21, 22, 23, 24, 25 }
	 M = 5; N = 5;
      case { 26, 27, 28, 29, 30 }
	 M = 5; N = 6;
      case { 31, 32, 33, 34, 35, 36 }
	 M = 6; N = 6;
      case { 37, 38, 39, 40, 41, 42 }
	 M = 6; N = 7;
      case { 43, 44, 45, 46, 47, 48, 49 }
	 M = 7; N = 7;
      case { 50, 51, 52, 53, 54, 55, 56 }
	 M = 7; N = 8;
      case 63
	 M = 9; N = 7;
      case { 57, 58, 59, 60, 61, 62, 64 }
	 M = 8; N = 8;
      otherwise
	 M = ceil(sqrt(nz)); N = ceil(sqrt(nz));
   end
end

return;
