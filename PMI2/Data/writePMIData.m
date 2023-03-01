% WRITEPMIDATA  Given SD and Data structures, write out a 
%               PMI-formatted data file
%
% writePMIData(SD, Data, filename);
%
% filename -> name of file to write.  The '.pmi' tag will be added
%              automatically

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

function writePMIData(SD, Data, filename)

% Force a .pmi extension, regardless of whether it was there originally
filename = strrep(filename, '.pmi', '');
filename = [ filename, '.pmi' ];

% Check that Data.PhiMeas and SD.Measlist match

if (~isfield(SD,'MeasList') | ~isfield(Data,'PhiMeas'))
  error('Incompatible input structures');
end

if (size(SD.MeasList,1) ~= size(Data.PhiMeas,1))
  error('Size mis-match between SD.MeasList and Data.PhiMeas');
end

% Open the output file, initially in text mode

fid = fopen(filename, 'wt');

if (fid < 0)
  error(['Error opening file ' filename]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write out fields as extracted from SD structure

fprintf(fid,'%% File generated from PMIDataFormat.m\n\n');

% Sources and detectors always get written out

for k = 1:size(SD.SrcPos,1)
  fprintf(fid, 'SrcPos(%3d) = [ %12.4e %12.4e %12.4e ]\n', k, ...
          SD.SrcPos(k,1), SD.SrcPos(k,2), SD.SrcPos(k,3));
end

fprintf(fid,'\n');

for k = 1:size(SD.DetPos,1)
  fprintf(fid, 'DetPos(%3d) = [ %12.4e %12.4e %12.4e ]\n', k, ...
          SD.DetPos(k,1), SD.DetPos(k,2), SD.DetPos(k,3));
end

fprintf(fid,'\n');

% Remainig fields are optional (write out only if defined)

if (isfield(SD, 'ModFreq')) % Frequency
  for k = 1:length(SD.ModFreq)
    fprintf(fid, 'Freq(k) = %f\n', SD.ModFreq(k));
  end
end

fprintf(fid,'\n');

if (isfield(SD, 'Lambda')) % ExcitationWavelength
  for k = 1:length(SD.Lambda)
    fprintf(fid, 'Lambda(k) = %f\n', SD.Lambda(k));
  end
end

fprintf(fid,'\n');

if (isfield(SD, 'EmissionWavelength'))
  for k = 1:length(SD.EmissionWavelength)
    fprintf(fid, 'EmissionWavelength(k) = %f\n', SD.EmissionWavelength(k));
  end
end

fprintf(fid,'\n');

if (isfield(SD, 'TimeDelay'))
  for k = 1:length(SD.TimeDelay)
    fprintf(fid, 'TimeDelay(k) = %f\n', SD.TimeDelay(k));
  end
end

fprintf(fid,'\n');

if (isfield(SD, 'TimeGateWidth'))
  for k = 1:length(SD.TimeGateWidth)
    fprintf(fid, 'TimeGateWidth(k) = %f\n', SD.TimeGateWidth(k));
  end
end

fprintf(fid,'\n');

if (isfield(SD, 'CorrelationTime'))
  for k = 1:length(SD.CorrelationTime)
    fprintf(fid, 'CorrelationTime(k) = %f\n', SD.CorrelationTime(k));
  end
end

fprintf(fid,'\n');

if (isfield(SD, 'DataType'))
  for k = 1:length(SD.DataType)
    fprintf(fid, 'DataType{k} = ''%s''\n', SD.DataType{k});
  end
end

fprintf(fid,'\n');

% etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out the measurement list

mlst = find(max(SD.MeasList) > 1);

for k = 1:size(SD.MeasList,1)
  fprintf(fid, 'Meas(%3d) = [ ', k);
  fprintf(fid, '%3d ', SD.MeasList(k, mlst));  % Vector write
  fprintf(fid, ']\n');
end

fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out a few final fields

fprintf(fid, 'DataPrecision = ''float64''\n');

if (~isfield(Data, 'DataType'))
   if (isfield(Data,'PhiMeas') & ~isreal(Data.PhiMeas))
      fprintf(fid, 'DataType(1)   = { ''IQ'' }\n');
   else			      
      fprintf(fid, 'DataType(1)   = { ''Amplitude'' }\n');
   end
end

fprintf(fid, '\n');
fprintf(fid, 'BeginData\n');

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the binary data

% Reopen the output file in binary mode

fid = fopen(filename, 'ab');

if (fid < 0)
  error(['Error re-opening file ' filename]);
end

count = fwrite(fid, Data.PhiMeas, 'float64');

fclose(fid);

return;


