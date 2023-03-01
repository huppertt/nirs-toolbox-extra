% CALCAFFINE  Find affine transform from (x,y) to (u,v),
%              i.e., find A,B such that [u v]' = A * [x y] + B
%
% [A,B] = calcAffine(x,y,u,v);

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

function[A,B] = calcAffine(x,y,u,v);

A = zeros(6,6);
B = zeros(6,1);

if (size(x) ~= size(y) | size(x) ~= size(u) | size(x) ~= size(v))
   error('Input vectors must all be the same size');
end

% Compute the individual matrix terms

sx  = sum(x);     sy  = sum(y);     su  = sum(u);     sv  = sum(v);
sxx = sum(x.*x);  sxy = sum(x.*y);  syy = sum(y.*y);  
sux = sum(u.*x);  suy = sum(u.*y);
svx = sum(v.*x);  svy = sum(v.*y);
s   = length(x);  % length(x) or 1 ???

% Construct the LLS matrix, form never changes

A(1,:) = [ sxx sxy   0   0 sx  0 ];
A(2,:) = [ sxy syy   0   0 sy  0 ];
A(3,:) = [   0   0 sxx sxy  0 sx ];
A(4,:) = [   0   0 sxy syy  0 sy ];
A(5,:) = [  sx  sy   0   0  s  0 ];
A(6,:) = [   0   0  sx  sy  0  s ];

% Construct the vector, again form is constant

B(:,1) = [ sux suy svx svy su sv ]';

% Find the least-squares solution

C = A \ B;

clear A B

A(1,:) = [ C(1) C(2) ];
A(2,:) = [ C(3) C(4) ];
B(:,1) = [ C(5) C(6) ]';

return;
