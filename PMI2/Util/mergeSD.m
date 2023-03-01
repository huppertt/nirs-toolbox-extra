% MERGESD Combine multiple SD structures into a single structure while
%         eliminating duplicate entries and preserving the ordering of
%         the measurement list.
%
% SD = mergeSD(SD1, SD2, ...);

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

function[SD] = mergeSD(varargin)

if (nargin == 1)
   % Remove duplicate entries and sort the fields of a single structure.
   
   SD = varargin{1};

   SD2 = SD;
   
   % Take out unsortable fields
   
   SD2.MeasList = [];
   
   if (isfield('DataType','SD2'))
      SD2 = rmfield(SD2,'DataType');
   end

   if (isfield('ImagerOption','SD2'))
      SD2 = rmfield(SD2,'ImagerOption');
   end
   
   % Pack down the structure
   
   SD = real_mergeSD(SD, SD2);
   
   return;
else
   SD = real_mergeSD(varargin{1}, varargin{2});
   
   for n = 3:nargin
      SD = real_mergeSD(SD, varargin{n});
   end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The real mergeSD routine, operates only on pairs of SD structures.

function[SD] = real_mergeSD(SD1, SD2);

% Get field names found in either structure

n1 = fieldnames(SD1);
n2 = fieldnames(SD2);

fn = unique([ n1; n2 ]);

map = cell(2,9);

ML1 = SD1.MeasList;
ML2 = SD2.MeasList;

% Initialize SD just to preserve legibility

SD.MeasList = [];
SD.SrcPos   = [];
SD.DetPos   = [];
SD.Lambda   = [];

% Merge the two SD structures

for k = 1:length(fn)
   fname = char(fn(k));
   
   if (isfield(SD1, fname) & ~isfield(SD2, fname))
      setfield(SD, fname, getfield(SD1, fname));
   elseif (isfield(SD2, fname) & ~isfield(SD1, fname))
      setfield(SD, fname, getfield(SD2, fname));
   elseif (isfield(SD2, fname) &  isfield(SD1, fname))
      if (strcmp(fname, 'MeasList'))
	 continue;          % Done on the fly
      elseif (strcmp(fname, 'SrcAmp') | strcmp(fname, 'DetAmp'))
	 continue;          % Deferred
      elseif (strcmp(fname, 'SrcPos'))
	 [M, ii, jj] = unique([ SD1.SrcPos; SD2.SrcPos ], 'rows');
	 map{1,1} = jj(1:size(SD1.SrcPos,1));
	 map{2,1} = jj(size(SD1.SrcPos,1)+1:end);
	    
	 [SD.SrcPos, ML1, ML2] = mapMatrix(ML1, ML2, 1, M, map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'DetPos'))
	 [M, ii, jj] = unique([ SD1.DetPos; SD2.DetPos ], 'rows');
	 map{1,2} = jj(1:size(SD1.DetPos,1));
	 map{2,2} = jj(size(SD1.DetPos,1)+1:end);
	 
	 [SD.DetPos, ML1, ML2] = mapMatrix(ML1, ML2, 2, M, map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'ModFreq'))
	 [M, ii, jj] = unique([ SD1.ModFreq(:)' SD2.ModFreq(:)' ]);
	 map{1,3} = jj(1:length(SD1.ModFreq));
	 map{2,3} = jj(length(SD1.ModFreq)+1:end);
	 
	 [SD.ModFreq, ML1, ML2] = mapMatrix(ML1, ML2, 3, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'Lambda'))
	 [M, ii, jj] = unique([ SD1.Lambda(:)' SD2.Lambda(:)' ]);
	 map{1,4} = jj(1:length(SD1.Lambda));
	 map{2,4} = jj(length(SD1.Lambda)+1:end);
	 
	 [SD.Lambda, ML1, ML2] = mapMatrix(ML1, ML2, 4, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'EmissionWavelength'))
	 [M, ii, jj] = unique([ SD1.EmissionWavelength(:)' ...
		                SD2.EmissionWavelength(:)' ]);
	 map{1,5} = jj(1:length(SD1.EmissionWavelength));
	 map{2,5} = jj(length(SD1.EmissionWavelength)+1:end);
	 
	 [SD.EmissionWavelength, ML1, ML2] = mapMatrix(ML1, ML2, 5, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'TimeDelay'))
	 [M, ii, jj] = unique([ SD1.TimeDelay(:)' SD2.TimeDelay(:)' ]);
	 map{1,6} = jj(1:length(SD1.TimeDelay));
	 map{2,6} = jj(length(SD1.TimeDelay)+1:end);
	 
	 [SD.TimeDelay, ML1, ML2] = mapMatrix(ML1, ML2, 6, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'TimeGateWidth'))
	 [M, ii, jj] = unique([ SD1.TimeGateWidth(:)'; ...
		                SD2.TimeGateWidth(:)' ]);
	 map{1,7} = jj(1:length(SD1.TimeGateWidth));
	 map{2,7} = jj(length(SD1.TimeGateWidth)+1:end);
	 
	 [SD.TimeGateWidth, ML1, ML2] = mapMatrix(ML1, ML2, 7, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'CorrelationTime'))
	 [M, ii, jj] = unique([ SD1.CorrelationTime(:)' ...
		                SD2.CorrelationTime(:)' ]);
	 map{1,8} = jj(1:length(SD1.CorrelationTime));
	 map{2,8} = jj(length(SD1.CorrelationTime)+1:end);
	 
	 [SD.CorrelationTime, ML1, ML2] = mapMatrix(ML1, ML2, 8, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'DataType'))
	 [M, ii, jj] = unique([ SD1.DataType(:)' SD2.DataType(:)' ]);
	 map{1,9} = jj(1:length(SD1.DataType));
	 map{2,9} = jj(length(SD1.DataType)+1:end);
	 
	 [SD.DataType, ML1, ML2] = mapMatrix(ML1, ML2, 9, M', map);
	 
	 clear M ii jj;
      elseif (strcmp(fname, 'ImagerOption'))
	 SD.ImagerOption = [ SD1.ImagerOption SD2.ImagerOption ];	 
      else
	 f1 = getfield(SD1, fname);
	 f2 = getfield(SD2, fname);
	 
	 % Not in measurement list, can simply append
	 setfield(SD, fname, [ f1(:)'; f2(:)' ]);
	 
	 clear f1 f2;
      end
   else
      error([ 'Impossible condition, field not found ' fname ]);
   end
   
   clear fname;
end

% Copy over the remapped measurement lists

SD.MeasList = [ ML1; ML2 ];

clear ML1 ML2;

% SD.SrcAmp and SD.DetAmp require some additional work

if (0)
SD.SrcAmp = ones(size(SD.SrcPos,1), ...
		 max(1,length(SD.Lambda)), max(1,length(SD.ModFreq)));

if (isfield(SD1,'SrcAmp') & isfield(SD1, 'MeasList'))
   for k = 1:size(SD1.MeasList,1)
      si = SD1.MeasList(k,1);
      fi = SD1.MeasList(k,3);
      wi = SD1.MeasList(k,4);
      
      SD.SrcAmp(map{1,1}(si),map{1,4}(wi),map{1,3}(fi)) = SD1.SrcAmp(si,wi,fi);
   end
end

if (isfield(SD2,'SrcAmp') & isfield(SD2, 'MeasList'))
   for k = 1:size(SD2.MeasList,1)
      si = SD2.MeasList(k,1);
      fi = SD2.MeasList(k,3);
      wi = SD2.MeasList(k,4);
      
      SD.SrcAmp(map{2,1}(si),map{2,4}(wi),map{2,3}(fi)) = SD2.SrcAmp(si,wi,fi);
   end
end

clear si fi wi k;

SD.DetAmp = ones(size(SD.DetPos,1), ...
		 max(1,length(SD.Lambda)), max(1,length(SD.ModFreq)));

if (isfield(SD1,'DetAmp') & isfield(SD1, 'MeasList'))
   for k = 1:size(SD1.MeasList,1)
      di = SD1.MeasList(k,2);
      fi = SD1.MeasList(k,3);
      wi = SD1.MeasList(k,4);
      
      SD.DetAmp(map{1,2}(di),map{1,4}(wi),map{1,3}(fi)) = SD1.DetAmp(di,wi,fi);
   end
end

if (isfield(SD2,'DetAmp') & isfield(SD2, 'MeasList'))
   for k = 1:size(SD2.MeasList,1)
      di = SD2.MeasList(k,2);
      fi = SD2.MeasList(k,3);
      wi = SD2.MeasList(k,4);
      
      SD.DetAmp(map{2,2}(di),map{2,4}(wi),map{2,3}(fi)) = SD2.DetAmp(di,wi,fi);
   end
end
end

clear di fi wi k;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map to new matrix.  This isn't a big function, but one copy is
% easier to debug than multiple "identical" copies.

function [M, ML1, ML2] = mapMatrix(ML1, ML2, nf, M, map);

tmp1 = map{1,nf}(ML1(:,nf));
tmp2 = map{2,nf}(ML2(:,nf));

ML1(:,nf) = tmp1(:);
ML2(:,nf) = tmp2(:);

return;

