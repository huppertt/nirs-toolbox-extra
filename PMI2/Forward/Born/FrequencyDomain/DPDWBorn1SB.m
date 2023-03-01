% DPDWBORN1SB    Forward DPDW model using 1st born for a slab medium.
%
%   One z boundary is assumed to be at z=0.  The other is assumed
%   to be at z = Model.Slab_Thickness
%
% [Phi0, A] = DPDWBorn1SB(SD, Medium, MeasList, Method, OptProp, Debug)
%
%   SD, Medium  The PMI Model structures
%   MeasList    The measurement list for this simulation.
%   Debug       OPTIONAL: Print out debugging info.
%
% Returns:
%   A           The forward matrix relating the contribution from each voxel to
%               a specific source-detector pair.  Each row is for a different
%               source detector pair.
%
%   Phi0        The incident response at each detector from each source
%               in a column vector.  The same combination pattern as the
%               columns of A.

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

function [Phi0, A] = DPDWBorn1SB(SD, Medium, MeasList, ...
                                   Method, OptProp, Debug)

if (~exist('MeasList','var') | isempty(MeasList))
   MeasList = SD.MeasList;
end

if (~exist('Debug','var'))
   Debug = 0;
end

%%
%%  Calculate the Born approximation and incident fields
%%

if strcmpi(Method,'Born') | strcmpi(Method,'Rytov')

   [Phi0, A] = Hlm3ptBorn1SB(SD, Medium, MeasList, OptProp, Debug);

   %%
   %%  If Method = 'Rytov' then normalize the matrix by Phi0
   %%

   if strcmp(Method,'Rytov')
      if (Debug)
        disp('Converting to Rytov form');
      end

      A = A ./ (Phi0 * ones(1,size(A,2)));
   end
elseif (strcmpi(Method,'FullBorn') | strcmpi(Method,'BornN')    | ...
	strcmpi(Method,'ExtBorn')  | strcmpi(Method,'ExtBornN'))

	warning('Badly implemented');

   [Phi0, A] = HlmFullBorn_slab2(SD, Medium, MeasList, OptProp, Debug);

%   % This part really belongs in HlmFullBorn_slab as well, but the
%   % incident field is harder to interpret in that case, so I'll leave
%   % it out here until I have a better idea how to implement it.
%   
%   wlst = MeasList(:,4);
%   
%   disp('DPDWBORN1SB: Hmmm... scaling for FullBorn right?');
%   Phi0 = -Phi0 * (Medium.v(1) / D(1));
%%   Phi0 = -Phi0 .* (Medium.v(wlst) ./ D(wlst))';
%   %   A = -(Medium.v(wlst) / D(wlst)) * A;
%   A = -A;
end

return;