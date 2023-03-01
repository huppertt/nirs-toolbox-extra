%function [mua] = GetMua( lambda, Cphore, Conc )
%
% Returns the absorption coefficient at the specified wavelengths
% given the chromophores and corresponding concentrations.
%
% Cphore is a structure array with fields:
%   Name - equals 'HbO', 'Hb', 'H2O', 'lipid', or 'aa3'
%   Conc - given in Molar except for water which is given as a
%   fraction.
%
% Conc is an optional vector and will contain the corresponding
%   concentrations rather than using those specified in the Cphore
%   structure array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: GetMua.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.1  1999/11/18 14:21:00  tgaudett
%  Initial Routines for Chromophores
%  HB, HbO, H2O
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mua] = GetMua( lambda, Cphore, Conc )

if nargin==2
   for i=1:length(Cphore)
      Conc(i) = Cphore(i).Conc;
   end
end

e(:,:) = feval('GetExtinctions', lambda );

num_lambda = length(lambda);
mua = zeros(size(Conc,1),num_lambda);

for i=1:length(Cphore)
   if strcmpi(Cphore(i).Name, 'HbO')
      mua = mua + Conc(:,i) * e(:,1)';
   elseif strcmpi(Cphore(i).Name, 'Hb')
      mua = mua + Conc(:,i) * e(:,2)';
   elseif strcmpi(Cphore(i).Name, 'H2O')
      mua = mua + Conc(:,i) * e(:,3)';
   elseif strcmpi(Cphore(i).Name, 'lipid')
      mua = mua + Conc(:,i) * e(:,4)';
   elseif strcmpi(Cphore(i).Name, 'aa3')
      mua = mua + Conc(:,i) * e(:,5)';
   end
end

%mua = mua';