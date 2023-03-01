%UVEC2STR       Convert a uniform vector to a string.
%
%   str = uvec2str(vec)
%
%   str         The resultant string.
%
%   vec         Teh uniform vector to be converted.
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
%  $Log: uvec2str.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%  Revision 1.2  1998/06/05 17:32:47  rjg
%  Handles 1 element vector now.
%
%  Revision 1.1  1998/04/28 20:36:20  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = uvec2str(vec)

if length(vec) > 1
    str = [ '[' num2str(vec(1)) ':' num2str(vec(2) - vec(1)) ':' ...
                num2str(vec(length(vec))) ']' ];
else
    str = [ '[' num2str(vec) ']' ];
end
