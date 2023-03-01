%cleanInvSys    Delete the forward system matrix from the PMI data structure.
%
%   pmi = cleanInvSys(pmi);
%
%   pmi          The PMI data structure to operate upon.
%
%
%   cleanInvSys removes inverse system to matrix to conserve memory.
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
%  $Log: cleanInvSys.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmi = cleanInvSys(pmi);

if isfield(pmi.Inv, 'A')
    pmi.Inv = rmfield(pmi.Inv, 'A');
end
if isfield(pmi.Inv, 'Aw')
    pmi.Inv = rmfield(pmi.Inv, 'Aw');
end


