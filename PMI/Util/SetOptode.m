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
%  $Date: 2000/07/27 15:13:38 $
%
%  $Revision: 1.2 $
%
%  $Log: SetOptode.m,v $
%  Revision 1.2  2000/07/27 15:13:38  dboas
%  The amplitude is now a 2 dimensional array.
%
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

function [pOptode] = SetOptode(X, Y, Z, Amp);

  [Xp Yp Zp] = meshgrid(X, Y, Z);
  pOptode.Pos = [Xp(:) Yp(:) Zp(:)];
  num = size(pOptode.Pos,1);
  
  pOptode.Type = 'list';
  if size(Amp,1) == num
    pOptode.Amplitude = Amp;
  else
    pOptode.Amplitude = ones(num,1) * Amp;
  end
