%calcInvModObjFct   Calculate the inverse model object function.
%
%   [delMu_a delMu_sp] = calcInvModObjFct(pmiStruct)
%
%   delMu_a     The absorption perturbation function evaluated at the inverse
%               computational model grid points.
%
%   delMu_sp    The scattering perturbation function evaluated at the inverse
%               computational model grid points.
%
%   pmiStruct   The PMI structure containing the Object structure and inverse
%               computational volume parameters.
%
%   calcInvModObjFct calculates the perturbation function(s) at each
%   wavelength returning the fuction evaluated at the inverse model
%   computation grid points.
%
%   Calls: calcDelMuA
%
%   Bugs: delMu_sp calculation not yet implemented

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: calcInvModObjFct.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.1  1999/12/22 16:04:05  rjg
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delMu_a] = calcInvModObjFct(pmi)

nLambda = length(pmi.Inv.Lambda);

for iLambda = 1:nLambda
    fct_cube = calcDelMuA(pmi.Object, pmi.Inv.CompVol, pmi.Inv.Mu_a, iLambda);

    delMu_a(:, iLambda) = fct_cube(:);
end
