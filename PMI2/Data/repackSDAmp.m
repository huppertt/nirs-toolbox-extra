%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the nSrc + nDet - 1 fitting coefficients, repack them into
%  vectors of source scaling and detector scaling coefficients.
%
% DO NOT CALL THIS FUNCTION DIRECTLY

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

function[svec, dvec] = repackSDAmp(fit, SD, MeasList)

nSrc = size(SD.SrcPos,1);
nDet = size(SD.DetPos,1);
nWav = length(SD.Lambda);

% Any term not in fit has a scale factor of 0.0 (Rytov)

svec = zeros(nSrc, nWav);
dvec = zeros(nDet, nWav);

offset = 0;

% Sources come first in the fit vector

for l = 1:nWav
   mlst = find(MeasList(:,4)==l);
   
   if ~isempty(mlst)
      vsrc = unique(MeasList(mlst,1));

      if (~isempty(vsrc))
	 svec(vsrc,l) = fit(offset + [1:length(vsrc)]);
	 offset = offset + length(vsrc); 
      end
      
      clear vsrc
   end
   
   clear mlst
end

% To reduce the rank of the matrix (which is initially
% rank-deficient), we deleted the first detector used at each
% wavelength, on a per-wavelength basis.  Figure out which
% coefficients this leaves me with and extract them from the fit
% vector.

for l = 1:nWav
   mlst = find(MeasList(:,4)==l);
   
   if ~isempty(mlst)
      vdet = unique(MeasList(mlst,2));
      vdet = vdet(2:end);    % First element removed to improve rank
      
      if (~isempty(vdet))
	 dvec(vdet,l) = fit(offset + [1:length(vdet)]);
	 offset = offset + length(vdet);
      end
      
      clear vdet
   end
   
   clear mlst 
end

return;

