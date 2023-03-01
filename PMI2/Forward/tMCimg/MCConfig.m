% MCCONFIG  Parse a computer friendly .inp (generated from a tMCimg
%           config file using cfg2inp)
%
% mc = MCConfig(filename);
%
% MCConfig() parses the given .inp file and returns a structure with the
%  different configuration parameters stored as structure fields.

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

function[mc] = MCConfig(filename)

fid=fopen([ filename '.inp' ], 'rt');
if (fid < 0)
   error([ 'Error opening file "' filename '.inp"' ]);
end

mc.NPh  = fscanf(fid,'%d', 1);

mc.seed = fscanf(fid,'%d', 1);

mc.freq = fscanf(fid,'%f', 1);
mc.freq = mc.freq * 1e6;       % Convert MHz to Hz

mc.minT   = fscanf(fid,'%f', 1);
mc.stepT  = fscanf(fid,'%f', 1);   
mc.nTstep = fscanf(fid,'%d', 1);
mc.maxT   = mc.minT + mc.stepT * mc.nTstep;

mc.binFile = fscanf(fid,'%s', 1);

mc.xstep = fscanf(fid,'%f', 1);
mc.nxstep = fscanf(fid,'%d', 1);
mc.ximin = fscanf(fid,'%f', 1);
mc.ximax = fscanf(fid,'%f', 1);

mc.ystep = fscanf(fid,'%f', 1);
mc.nystep = fscanf(fid,'%d', 1);
mc.yimin = fscanf(fid,'%f', 1);
mc.yimax = fscanf(fid,'%f', 1);

mc.zstep = fscanf(fid,'%f', 1);
mc.nzstep = fscanf(fid,'%d', 1);
mc.zimin = fscanf(fid,'%f', 1);
mc.zimax = fscanf(fid,'%f', 1);

mc.nTissue = fscanf(fid,'%d', 1);
for idx=1:mc.nTissue
   % mus, g, mua, n
   mc.tis(idx,:) = fscanf(fid,'%f', 4)';
end

% number, radius, direction - sources
mc.nSrc(1) = fscanf( fid, '%d', 1);
mc.nSrc(2) = fscanf( fid, '%f', 1);
mc.ci(:,1) = fscanf( fid, '%f', 3);

mc.sr  = mc.nSrc(2);
mc.sna = 0;

for idx=1:mc.nSrc(1)
   mc.xi(:,idx) = fscanf( fid, '%f', 3);
end

% number, radius - detector
mc.nDet(1) = fscanf( fid, '%d', 1);
mc.nDet(2) = fscanf( fid, '%f', 1);

for idx=1:mc.nDet(1)
   mc.pOpt(:,idx) = fscanf( fid, '%f', 3);
end

fclose(fid);

return;

