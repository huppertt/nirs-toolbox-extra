% DPDWHELMHOLTZ  Forward DPDW model using exact solution of Helmholtz equation
%
% [Phi0] = DPDWHelmholtz(SD, Medium, MeasList);
%
%   SD, Medium - PMI SD data structure
%   MeasList   - The measurement list for this simulation.
%
% Returns:
%   Phi0     The incident response at each detector from each source
%               in a column vector.  
%
% NOTE: DPDWHelmholtz() tries to guess whether to use a Frequency-Domain or
%   a Time-Domain Green's function based on which fields are defined in SD.
%   Call FD2pt() or TD2pt() explicitly if you get errors.  CW imaging is
%   treated as a special case of FD imaging.

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

function [Phi0] = DPDWHelmholtz(SD, Medium, MeasList, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var') | isempty(Debug))
   Debug = 0;
end

%% Decide what style imager this is based on the existance and value
%%  of certain fields in the SD structure and the Measurement list.  
%%  Fields whose value is an empty matrix ([]) are treated as missing.

is_td = 0;    % Initialize the flags to zero
is_rf = 0;

is_td = isTD(SD, MeasList);
is_rf = isFD(SD, MeasList);

if (is_rf == 0 & is_td == 0)
   % Give up now
   error('Unable to guess imager design');
end

%%
%%  Select the appropriate matrix generation routine
%%

if (is_td & is_rf)
   % This requires human intervention
   error(['Imager would appear to be both CW/RF _and_ TD!' ...
	  '  Call FD2pt or TD2pt directly.']);
elseif (is_td)
   % Call high-level TD wrapper
   Phi0 = TD2pt(SD, Medium, MeasList, Debug);
elseif (is_rf)
   % Call high-level FD wrapper
   Phi0 = FD2pt(SD, Medium, MeasList, Debug);
else
   % Not enough information to proceed
   error('Imager is neither CW, RF, nor TD!');
end

return;
