% GENMEASLIST    Generate a measurement list 
%
%   MeasList = genMeasList(SD, method, var1, var2, ...);
%
%   method      This defines the method used to create a subset
%               measurement list.  Presently it can equal
%                 'Nearest' - Choose the nearest neighbor source-detector
%                 	  pairs such that the separation is less than var1.
%	 	  'SNR' - Choose the threshold SNR.
%			  var1 - threshold SNR
%			  var2 - data frame
%                         var3 - std. deviation of data frame
%                 'All' - Choose all combinations of source,
%                         detectors, modulation frequencies, and
%                         wavelengths.
%                 'Range' Choose the neighbor source-detector pairs such
%                         that the separation is within given range
%                         between var1 and var2
%
%   var1, var2, etc
%               Variables whose use depends on METHOD.
%
% Returns:
%   MeasList - The measurement list in the PMI datastructure format.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2004, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang,
%                     Jonathan Stott, Daniel Haensse
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

function [MeasList] = genMeasList(SD, method, varargin)

% Do not use EmissionWavlength (really SD.Lambda), CorrelationTime (really
% SD.TimeDelay) or DataType in automatically generated measurement lists.
% These fields will probably require a custom Method at some point in the
% future.

MLField = { 'SrcPos', 'DetPos', 'ModFreq', 'Lambda', ...
	    'XEmissionWavelength', 'TimeDelay', 'TimeGateWidth', ...
	    'XCorrelationTime', 'XDataType' };

% SrcPos and DetPos are required fields, all others are optional

if (~isfield(SD,'SrcPos') | isempty(SD.SrcPos) | ...
    ~isfield(SD,'DetPos') | isempty(SD.DetPos))
   error('Both SD.SrcPos and SD.DetPos must exist');
else
   ne(1) = size(SD.SrcPos,1);
   ne(2) = size(SD.DetPos,1);
end

nField = length(MLField);

for k = 3:nField
   if (isfield(SD, char(MLField(k))))
      ne(k) = length(getfield(SD,char(MLField(k))));
   else
      ne(k) = 1;
   end
end

nMeas = prod(ne);

MeasList = zeros(nMeas, nField);

% Set up a list with all possible combinations as a starting point

stride = 1;

for k = 1:nField
   li = 1;
   
   if (isfield(SD, char(MLField(k))))
      for l = 1:stride:nMeas
	 i1 = l;
	 i2 = l + stride - 1;
	 
	 MeasList(i1:i2, k) = mod(li-1, ne(k)) + 1;
	 li = li + 1;
      end
   end
   
   stride = stride * ne(k);
end

% Throw out undesirable points

switch lower(method)
   case 'range'
      R1 = double(varargin{1});
      R2 = double(varargin{2});
      dR = calcSep(SD, MeasList);

      ml = find (dR > R1 & dR <= R2)

      MeasList = MeasList(ml, :);

   case 'nearest'
      R  = double(varargin{1});
      dR = calcSep(SD, MeasList);
      
      ml = find(r <= R);
  
      MeasList = MeasList(ml, :);
  
   case 'snr'
      SNRx = double(varargin{1});
      data = double(varargin{2});
      std  = double(varargin{3});

      snr = abs(data ./ std);
      ml  = find(snr > SNRx);

      MeasList = MeasList(ml,:);
   
   case 'all'
      MeasList = MeasList(1:end,:);
   
   otherwise
      warning([ 'Unknown method "' method ...
		'" for generating a measurement list.']);
end

return;
