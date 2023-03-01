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

function [phi_inc, phi_sc] = ExactSphere( rs, robj, rd, ko, Do, ki, ...
					  Di, a, vo)
  rso = sum((rs-robj).^2)^0.5;
  rod = sum((rd-robj).^2)^0.5;
  rsd = sum((rs-rd).^2)^0.5;
  c_theta = ((robj-rs) * (rd-robj)') / (rso*rod);
    
  x = ko*a;
  y = ki*a;
  
  phi_sc = 0;

  for l=0:20

    SHs = SphHarm(l,-1,0);
    SHd = SphHarm(l,c_theta,0);
    
    A = -j * ko * bes_h1(l,ko*rso) * ...
	( (Do*x*bes_jd(l,x)*bes_j(l,y) - Di*y*bes_j(l,x)*bes_jd(l,y)) ...
	 /(Do*x*bes_h1d(l,x)*bes_j(l,y)- Di*y*bes_h1(l,x)*bes_jd(l,y)) ...
	  );

    foo = SHs(1) * SHd(1);
    for m=1:l
      foo = foo + conj(SHs(m+1)) * SHd(m+1);
      foo = foo + SHs(m+1) * conj(SHd(m+1));
    end
    
    foo1 = A * bes_h1(l,ko*rod);
    phi_sc = phi_sc + foo1 * foo;

  end
  
%  phi_3pt = (vo/Do) * (exp(j*ko*rso+j*ko*rod)/(4*pi*rso*4*pi*rod)) * (4/3*pi*a^3*vo/Do);
  phi_inc = vo * exp(j*ko*rsd)/(4*pi*Do*rsd);
  phi_sc = vo/Do * phi_sc;
    