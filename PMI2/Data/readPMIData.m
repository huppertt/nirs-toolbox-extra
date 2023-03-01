% READPMIDATA   read in and parse a PMI structured datafile
%
% [SD, Data] = readPMIData(filename);
%
% Inputs: filename - name of the PMI structured datafile
%
% Outputs: SD   - PMI data structures
%          Data - data table read from datafile

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

function[SD, Data] = readPMIData(filename)

Debug = 1;

% Initialize the core fields in SD structure

SD.Lambda   = [];
SD.SrcPos   = [];
SD.DetPos   = [];
SD.MeasList = [];
SD.DataPrec = 'real*8';

% Open the file for input

fid = fopen(filename,'rb');

if (fid < 0)
   warning(['Error opening file ' filename]);
   keyboard;
end

SD = parseHeader(SD, fid, filename);

% Read the data portion of the file ------------------

ndata  = size(SD.MeasList,1);
nframe = 1;
Data   = [];
cflag  = 0;

if (isfield(SD, 'DataType'))
   if (length(SD.DataType) > 1)
      warning('Multiple SD.DataTypes are not yet supported');
      keyboard;
   end
   
   if (strcmpi(SD.DataType{1}, 'IQ') | ...
       strcmpi(SD.DataType{1}, 'Complex'))

      % Now it gets fun.  Matlab won't let me read or write complex
      % data types directly, but CW/RF filters generate IQ
      % (i.e. complex) data, so we need some way to handle this.  What
      % we do now is set the complex flag, double the size of our
      % reads, and repack as complex before copying back each frame of
      % data.  Someday this will need to be extended to handle mixed
      % DataType's, but not right now.
      
      cflag = 1;
   end
end

% No reason to read data portion if we're not going to return it to the
% user.

if (nargout > 1)
   % Preserve data type across reads.  This is MUCH more efficient than
   % converting everything into a double when dealing with large file
   % sizes.

   [frames, count] = fread(fid, inf, [ SD.DataPrec ' => ' SD.DataPrec ]);
   
   if (mod(count, ndata) ~= 0)
      warning('Read incomplete frames--truncated file?');
      keyboard;
   
      % In case the user insists on continuing
      if (cflag)
	 count = count - mod(count, 2*ndata);
      else
	 count = count - mod(count,   ndata);
      end
      
      frames = frames(1:count);
   end
   
   if (cflag)
      if (mod(count,2) ~= 0)
	 warning('Complex data, but odd number of floats read');
	 keyboard;
      end
   
      % Make frames complex and adjust count.  Only double-complex data
      % types are supported by Matlab, so we have no choice but to promote
      % here.  I still think it's faster, though, to read in as single and
      % promote afterwards than to let fread do the conversion.
   
      frames = double(frames(1:2:count)) + i*double(frames(2:2:count));
      count = count / 2;
   end

   % Pack the data frame into the PMI structure

   Data = reshape(frames, ndata, count / ndata);
   nframes = size(Data,2);
   clear frames;
end

% Done with file
fclose(fid);

% DataPrec is no longer useful, the data has either been promoted to double
%  by matlab or DataPrec can be discovered using just "whos".

SD = rmfield(SD, 'DataPrec');

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the header portion of the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[SD] = parseHeader(SD0, fid, filename)

SD = SD0;

MeasList = [];
txtline  = 'xyzzy';		% Get into the while() loop

while (~strcmp(txtline,'BeginData'))
   txtline = fgetl(fid);
   cmdline = txtline;
   
   % Strip off comments
   
   ci = findstr(cmdline, '%');
   
   if (~isempty(ci))
      cmdline = cmdline(1:ci(1)-1);
   end
   
   % Parse this line of text

   keyword = strtok(cmdline, ' =(){}');

   if     (isempty(keyword))
      % blank line, do nothing
   elseif (strcmp(keyword, 'BeginData'))
      % Switch to reading binary data
      break;
   else
      % Keyword processing
      
      switch keyword
	 % MeasurementList fields
	 
	 case { 'SrcPos' }
	    cmdline = strrep(cmdline, keyword, 'SD.SrcPos');
	    cmdline = strrep(cmdline, ')', ',:)');
	    eval([ cmdline ';' ],'parseError(filename, txtline)');
	 case { 'DetPos' }
	    cmdline = strrep(cmdline, keyword, 'SD.DetPos');
	    cmdline = strrep(cmdline, ')', ',:)');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	 case { 'Frequency' 'ModFreq' }
	    cmdline = strrep(cmdline, keyword, 'SD.ModFreq');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	 case { 'Lambda' 'ExcitationWavelength' }
	    cmdline = strrep(cmdline, keyword, 'SD.Lambda');
	    eval([ cmdline ';' ], 'parseError(filename, cmdline)');
	 case { 'EmissionWavelength' }
	    cmdline = strrep(cmdline, keyword, 'SD.EmissionWavelength');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	 case { 'TimeDelay' }
	    cmdline = strrep(cmdline, keyword, 'SD.TimeDelay');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	 case { 'TimeGateWidth' }
	    cmdline = strrep(cmdline, keyword, 'SD.TimeGateWidth');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	 case { 'CorrelationTime' }
	    cmdline = strrep(cmdline, keyword, 'SD.CorrelationTime');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	    
	    % Items in the measurement list
	    
	 case { 'Meas' }
	    cmdline = strrep(cmdline, keyword, 'MeasList');
	    cmdline = strrep(cmdline, ')', ',:)');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	    
	    % Data formatting options
	    
	 case { 'DataPrecision' }
	    cmdline = strrep(cmdline, keyword, 'SD.DataPrec');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	 case { 'DataType' }
	    if (isfield(SD,'DataType'))
	       warning('SD.Datatype multiply defined -- not supported');
	       keyboard;
	    end
	    
	    cmdline = strrep(cmdline, keyword, 'SD.DataType');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	    
	    % Misc fields
	    
	 case { 'ImagerOption' }
	    cmdline = strrep(cmdline, keyword, 'SD.ImagerOption');
	    eval([ cmdline ';' ], 'parseError(filename, txtline)');
	    
	    % Error - unknow fields
	    
	 otherwise
	    warning([ 'Unknown keyword "' keyword '"' ]);
	    keyboard;
      end
   end
end

SD.MeasList = repackMeasList(SD, MeasList);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repack the measurements into a proper measurement list
%
% Missing fields get filled in with zeros.
%
% Fields with a single value are implicitly 1, but were
%  were not specified in the datafile
%
% Fields with a vector value were explicitly specified and
%  can be copied directly into the appropriate column.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[MeasList] = repackMeasList(SD, ML)

nfields = 9;              % Width of measurement list
mi = 1;                   % index into datafile measurement columns

% Initialize blank fields to all zeros

SD.MeasList = zeros(size(ML,1), nfields);

% Pack the specified columns into the measurement list

if (isfield(SD,'SrcPos') & ~isempty(SD.SrcPos))
   if (size(SD.SrcPos,1) == 1)              % Field 1, source index
      SD.MeasList(:,1) = 1;
   else
      SD.MeasList(:,1) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'DetPos') & ~isempty(SD.DetPos))
   if (size(SD.DetPos,1) == 1)               % Field 2, detector index
      SD.MeasList(:,2) = 1;
   else
      SD.MeasList(:,2) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'ModFreq') & ~isempty(SD.ModFreq))
   if (length(SD.ModFreq) == 1)               % Field 3, frequency
      SD.MeasList(:,3) = 1;
   else
      SD.MeasList(:,3) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'ExcitationWavelength') & ~isempty(SD.ExcitationWavelength))
   if (length(SD.ExcitationWavelength) == 1) % Field 4a, excitation
      SD.MeasList(:,4) = 1;
   else
      SD.MeasList(:,4) = ML(:,mi);
      mi = mi + 1;
   end
elseif (isfield(SD, 'Lambda') & ~isempty(SD.Lambda))
   if (length(SD.Lambda) == 1)               % Field 4b, excitation
      SD.MeasList(:,4) = 1;
   else
      SD.MeasList(:,4) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'EmissionWavelength') & ~isempty(SD.EmissionWavelength))
   if (length(SD.EmissionWavelength) == 1)   % Field 5, emission
      SD.MeasList(:,5) = 1;
   else
      SD.MeasList(:,5) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'TimeDelay') & ~isempty(SD.TimeDelay))
   if (length(SD.TimeDelay) == 1)            % Field 6, delay time
      SD.MeasList(:,6) = 1;
   else
      SD.MeasList(:,6) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'TimeGateWidth') & ~isempty(SD.TimeGateWidth))
   if (length(SD.TimeGateWidth) == 1)        % Field 7, gate width
      SD.MeasList(:,7) = 1;
   else
      SD.MeasList(:,7) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'CorrelationTime') & ~isempty(SD.CorrelationTime))
   if (length(SD.CorrelationTime) == 1)      % Field 8, correlation time
      SD.MeasList(:,8) = 1;
   else
      SD.MeasList(:,8) = ML(:,mi);
      mi = mi + 1;
   end
end

if (isfield(SD, 'DataType') & ~isempty(SD.DataType))
   if (length(SD.DataType) == 1)             % Field 9, data type
      SD.MeasList(:,9) = 1;
   else
      SD.MeasList(:,9) = ML(:,mi);
      mi = mi + 1;
   end
end

if (mi ~= size(ML, 2)+1)
   warning('Fields defined but not found in measurement list');
   keyboard;
end

MeasList = SD.MeasList;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Error occured evaluating line, print message and stop.

function parseError(filename, txt)

error(['Error parsing file ' filename ', last line read was ' txt ]);

return;
