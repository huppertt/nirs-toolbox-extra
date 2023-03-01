%MAT2STR        Convert a vector to a string.
%
%   str = mat2str(vec)
%
%   str         The resultant string.
%
%   vec         The vector to be converted.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = mat2str(vec)

nRows = size(vec,1);
nCols = size(vec,2);
if nElem == 1
    str = num2str(vec);
    return;
end
if nRows == 1
    str = vec2str(vec);
end

str = '[';
for iRow = 1:Rows
    for iCol = 1:nCols
        str = [ str num2str(vec(iRow, iCol)) ' '];
    end
    str = [ str ';' ];
end
str = [str ']' ];

