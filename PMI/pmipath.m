%pmipath        Add the necessary PMI directories to the MATLAB path.
%
%   pmipath('Basepath');
%
%   Make sure to set the BasePath variable in this function to the base path at
%   your site.

% Copyright (C) 2002, David Boas, Dana Brooks, Rick Gaudette, 
%                     Tom Gaudette, Eric Miller, Quan Zhang
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: tgaudett $
%
%  $Date: 2000/05/25 14:27:52 $
%
%  $Revision: 1.2 $
%
%  $Log: pmipath.m,v $
%  Revision 1.2  2000/05/25 14:27:52  tgaudett
%
%  This needed to be done allong time ago.  All you need to pass the pmipath now is the base directory of where you put the repository.  So normally you go to the directory in matlab and type pmipath(pwd);
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.5  1999/11/18 17:41:14  dboas
%  Added instruments directory to path
%
%  Revision 3.4  1999/11/15 19:34:15  rjg
%  Added TimeSeries to the path.
%
%  Revision 3.3  1999/11/09 22:22:04  rjg
%  Removed non core directories
%
%  Revision 3.2  1999/10/21 21:56:44  rjg
%  Removed DOS linefeeds.
%
%  Revision 3.1  1999/10/19 21:14:02  rjg
%  Fixed DOS linefeeds.
%
%  Revision 3.0  1999/06/17 19:29:38  rjg
%  Initial Revision for PMI 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmipath(BasePath);

%%
%%  PMI directories
%%
if strcmp(computer, 'PCWIN')
    PMIPath = [ ...
            BasePath ';' ...
            BasePath '\Forward;' ...
            BasePath '\Reconstruct;' ...
            BasePath '\Reconstruct\PerfMeas;' ...
            BasePath '\TimeSeries;' ...
            BasePath '\Visualize;' ...
            BasePath '\Instruments;' ...
            BasePath '\UI;' ...
            BasePath '\Util' ];
else
    PMIPath = [ ...
            BasePath ':' ...
            BasePath '/Forward:' ...
            BasePath '/Reconstruct:' ...
            BasePath '/Reconstruct/PerfMeas:' ...
            BasePath '/TimeSeries:' ...
            BasePath '/Visualize:' ...
            BasePath '/Instruments:' ...
            BasePath '/UI:' ...
            BasePath '/Util' ];
end
            
path(path, PMIPath);


disp('===================================================================')
disp('PHOTON MIGRATION IMAGING TOOLBOX')
disp('Copyright (C) 2002, David Boas, Dana Brooks, Rick Gaudette,') 
disp('                    Tom Gaudette, Eric Miller, Quan Zhang')
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
disp('===================================================================')
