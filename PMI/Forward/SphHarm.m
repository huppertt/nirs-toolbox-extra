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

function foo = SphHarm( l, c_theta, phi )
  
  a = legendre(l,c_theta);
  for m=0:l
    foo(m+1) = sqrt( (2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m) ) * a(m+1) * exp(j*m*phi);
  end
  