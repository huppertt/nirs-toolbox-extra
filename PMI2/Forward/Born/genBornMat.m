% GENBORNMAT  Generate the Born Matrix for the specified method.  
%
% [Phi0, A] = genBornMat(SD, Medium, MeasList, Method, OptProp);
%   
%   Method - This specifies the method used to solve the Forward
%            Problem.  It can be one of the following:
%               'Born'  - Solves the first Born approximation,
%               'Rytov' - Solves the Rytov approximation,
%               'ModifiedBorn' - Time-domain analog of 'Rytov',
%               'FullBorn', 'BornN', 'ExtBorn', 'ExtBornN'
%                       - Higher order Born approximations.
%            Not all methods are supported by all imagers.
%
%   SD, Medium - PMI data structures
%   MeasList   - The Measurement List
%
%   OptProp - is a flag with 2 elements. 1-Yes, 0-No.
%      OptProp(1) indicates whether to consider perturbations in
%         the absorption coefficient when solving the Forward
%         Problem.
%      OptProp(2) indicates whether to consider perturbations in
%         the scattering coefficient when solving the Forward
%         Problem.
%
% Returns:
%    Phi0 - The total measured fluence for each
%         source-detector pair specified in SD.MeasList given a
%         homogneeous medium.  
%    A - The matrix used to calculate the Forward Problem.

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

function[Phi0, A] = genBornMat(SD, Medium, MeasList, Method, OptProp, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var'))
   Debug = 0;
end

%% Decide what style imager this is based on the existance and value
%%  of certain fields in the SD structure and the Measurement list.  
%%  Fields whose value is an empty matrix ([]) are treated as missing.

is_td = 0;    % Initialize the flags to zero
is_fd = 0;

if (isTD(SD,MeasList))
  is_td = 1;
end

if (isFD(SD,MeasList))
  % CW is a special case of frequency-domain
  is_fd = 1;
end

if (is_fd == 0 & is_td == 0)
   % Give up now
   error('Unable to guess imager design');
end

if (is_fd & is_td)
   error('Cannot handle mixed TD-FD designs');
end

%%
%%  Select the appropriate matrix generation routine
%%

if (is_fd)            % CW is just a variation on the RF imagers
   if (strcmpi(Method,'Born') | strcmpi(Method,'Rytov'))
      [Phi0, A] = FD3pt(SD, Medium, MeasList, OptProp, Debug);
      
      % Turn Born forward matrix into Rytov as needed

      if (strcmpi(Method,'Rytov'))
	 if (Debug)
	    disp('Converting from Born to Rytov form');
	 end

	 for k = 1:size(A,2)
	    A(:,k) = A(:,k) ./ Phi0;
	 end
      end
   elseif (strcmpi(Method,'FullBorn') | strcmpi(Method,'BornN') | ...
	   strcmpi(Method,'ExtBorn')  | strcmpi(Method,'ExtBornN'))
      
      switch lower(Medium.Geometry)
	 case {'infinite', 'inf', 'inft' }
	    [Phi0,A] = HlmFullBornNB(SD,Medium,MeasList,OptProp,Debug);
      
	 case { 'semi-infinite', 'semi' }
	    [Phi0,A] = HlmFullBornZB(SD,Medium,MeasList,OptProp,Debug);
      
	 case {'slab'}
	    [Phi0,A] = HlmFullBornSB(SD,Medium,MeasList,OptProp,Debug);
	    
	 otherwise
      end
   else
      error(['Unsupported FD method ' Method]);
   end
elseif (is_td)
   % The TD code will figure out the geometry, don't bother here
   
   if (strcmpi(Method,'Born') | strcmpi(Method,'ModifiedBorn'))
      % Currently, only Born method is supported
      [Phi0, A] = TD3pt(SD, Medium, MeasList, OptProp, Debug);
	    
      % Rytov doesn't carry over to time-domain, but amplitude-
      % weighted reconstructions are still useful
      
      if (strcmpi(Method,'ModifiedBorn'))
	 for k = 1:size(A,2)
	    A(:,k) = A(:,k) ./ Phi0;
	 end
      end
   elseif (strcmpi(Method,'FullBorn'))
      error([ 'TD does not yet support FullBorn' ]);
   else
      error(['Unknown TD method ' Method ]);
   end
else
   error('Unable to proceed - SD would appear to be neither TD nor FD');
end

return;
