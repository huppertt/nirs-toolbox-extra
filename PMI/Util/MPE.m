% MPE - This function calculates the maximum permissible tissue
%       exposure given the wavelength of light and the exposure
%       duration.  Above 10 s the MPE only depends on
%       wavelength. Below 10 s it has a t^(0.25) dependence.
%
% FluenceRate = MPE( lambda, time )
%
% lambda - The wavelength in nanometers.
% time   - The exposure time in seconds.
%
% FUNCTION RETURNS
% FluenceRate - This is the MPE expressed as a power density [Watts/mm^2]
% allowed during the given time interval.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/09/06 12:23:32 $
%
%  $Revision: 1.1 $
%
%  $Log: MPE.m,v $
%  Revision 1.1  2000/09/06 12:23:32  dboas
%  initial revision
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FluenceRate = MPE( lambda, time)
  
  Ca = 0.04*exp(4.597*lambda/1e3);
  
  if time < 10
    FluenceRate = 1.1 * Ca .* time.^0.25 ./ 100 ./ time;
  else
    FluenceRate = 0.2 * Ca ./ 100;
  end
  
  