%getOptodePos   Get the optode position and number from the optode data structure.
%
%   [pOptode, nOptode] = getOptodePos(OptodePos);
%
%   pOptode     A array of the positions of each optode.
%
%   noptode     The number of optodes present.
%
%   OptodePos   The PMI Optode data structure.
%
%
%   Calls: none.
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: getOptodePos.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.2  1999/11/10 19:05:26  rjg
%  Removed DOS linefeeds.
%
%  Revision 3.1  1999/10/13 17:38:24  rjg
%  Fixed pOptode assignment when the type is a list.
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pOptode, nOptode] = getOptodePos(OptodePos);
if strcmp(OptodePos.Type, 'uniform')
    [Xp Yp Zp] = meshgrid(OptodePos.X, OptodePos.Y,OptodePos.Z);
    nOptode = prod(size(Xp));
    pOptode = [Xp(:) Yp(:) Zp(:)];
else
    nOptode = size(OptodePos.Pos, 1);
    pOptode = OptodePos.Pos;
end
