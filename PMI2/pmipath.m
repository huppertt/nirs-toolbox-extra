% PMIPATH  Add the PMI toolbox directories to the MATLAB search path.
%
% pmipath('Basepath');
%
% Make sure to set the BasePath variable in this function to the
% base path at your site.  This progam must be run before using the
% rest of the toolbox (we recommend putting it directly in your Matlab
% startup.m file).

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

function pmipath(BasePath);

if (~exist('BasePath','var') | isempty(BasePath))
   error('No basepath specified')
end

%%
%% PMI directories, add at head of path (over-rides existing functions)
%% Using fullfile() makes this code OS independant.
%% 

addpath(BasePath, ...
	fullfile(BasePath, 'Data'), ...
	fullfile(BasePath, 'Forward'), ...
	fullfile(BasePath, 'Forward', 'Born', 'MEX'), ...
	fullfile(BasePath, 'Forward', 'Born'), ...
	fullfile(BasePath, 'Forward', 'Born', 'Correlation'), ...
	fullfile(BasePath, 'Forward', 'Born', 'Fluorescence'), ...
	fullfile(BasePath, 'Forward', 'Born', 'FrequencyDomain', 'MEX'), ...
	fullfile(BasePath, 'Forward', 'Born', 'FrequencyDomain'), ...
	fullfile(BasePath, 'Forward', 'Born', 'TimeDomain', 'MEX'), ...
	fullfile(BasePath, 'Forward', 'Born', 'TimeDomain'), ...
	fullfile(BasePath, 'Forward', 'Sphere'), ...
	fullfile(BasePath, 'Forward', 'tMCimg'), ...
	fullfile(BasePath, 'Forward', 'tFDimg'), ...
	fullfile(BasePath, 'Forward', 'TOAST'), ...
	fullfile(BasePath, 'Noise'), ...
	fullfile(BasePath, 'Recon'), ...
	fullfile(BasePath, 'Util'), ...
	fullfile(BasePath, 'Visualize'), ...
	fullfile(BasePath, 'Local'), ...
	fullfile(BasePath, 'Development'), '-begin' );
            
if (1)
   disp('===================================================================')
   disp('PHOTON MIGRATION IMAGING TOOLBOX')
   disp('Copyright (C) 2004, David Boas, Dana Brooks, Rick Gaudette, ')
   disp('                    Tom Gaudette, Eric Miller, Quan Zhang,')
   disp('                    Jonathan Stott')
   
   if (0)
      disp(' ')
      disp('This program is free software; you can redistribute it and/or')
      disp('modify it under the terms of the GNU General Public License')
      disp('as published by the Free Software Foundation; either version 2')
      disp('of the License, or (at your option) any later version.')
      disp(' ')
      disp('This program is distributed in the hope that it will be useful,')
      disp('but WITHOUT ANY WARRANTY; without even the implied warranty of')
      disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the')
      disp('GNU General Public License for more details.')
      disp(' ')
      disp('You should have received a copy of the GNU General Public License')
      disp('along with this program; if not, write to the Free Software')
      disp('Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA ')
      disp('                                          02111-1307, USA.')
   end
   
   disp('===================================================================')
end
