function [x,varargout] = GetConcentrations( cphores, lambda, muas )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: tgaudett $
%
%  $Date: 2000/08/17 15:58:15 $
%
%  $Revision: 1.2 $
%
%  $Log: GetConcentrations.m,v $
%  Revision 1.2  2000/08/17 15:58:15  tgaudett
%  Added the Ability to return residual if you want it.
%
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

num_lambda = length(lambda);
num_cphores = length(cphores);

ee(:,:) = GetExtinctions( lambda );

A = zeros(num_lambda, num_cphores);
for i=1:num_cphores;
   if strcmp(cphores(i).Name,'HbO')  
      A(:,i) = ee(:,1); end
   if strcmp(cphores(i).Name,'Hb')  
      A(:,i) = ee(:,2); end
   if strcmp(cphores(i).Name,'H2O')  
      A(:,i) = ee(:,3); end
end

% Solve the problem of b = Ax
%  inv(A)*b = x
%

if num_lambda==num_cphores
	% Then least squares does not buy you anything so just solve
	c=inv(A)*muas';
	if nargout>1
		disp('Error No Residual');
	end;	
	x= c';
else
	x=inv(A'*A)*A'*muas';
	resid_t = muas'-A*x;
	resnorm_t = sum(resid_t.*resid_t);
	varargout(1) = {resid_t};
	varargout(2) = {resnorm_t}; %This is the 2-norm of the Residual;
	x = x';
end;
return;
%
%  I needed to use names like resid_t instead of residual because people 
%  already used these names as function names and matlab has a resid 
%  function already.


