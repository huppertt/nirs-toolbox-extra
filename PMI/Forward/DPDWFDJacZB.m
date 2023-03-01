%DPDWFDJacZB    DPDW Half Space Finite Difference Jacobian Iteration Solution
%
%   [PhiScat PhiInc] = DPDWFDJacZB(CompVol, mu_sp, mu_a, delMu_a, v, f,
%                                  pSrc, Asrc, pDet, ADet, Debug)
%
%
%   PhiScat     The scattered field at each detector.  This vector is structured
%               as a column vector with detectors being iterated over first
%               followed by the sources.
%
%   PhiInc      The incident response at each detector from each source in a
%               column vector.
%
%   CompVol     A structure defining the computational volume.  This
%               structure should have the members: Type, X, Y and Z.  Type
%               should be uniform specifying a uniform sampling volume of
%               voxels. X, Y and Z are vectors specifying the centers of the
%               voxels.
%
%   mu_sp       The reduced scattering coefficient.
%
%   mu_a        The background absorption coeffiecient.
%
%   v           The propagation velocity in the medium.
%
%   f           The modulation frequency of the source.
%
%   pSrc        The position of the source(s) in the form [sx sy sz].
%               Each row represents a different source.
%
%   ASrc        The amplitude of each source.
%
%   pDet        The position of the detector(s) in the form [sx sy sz].
%               Each row represents a different detector.
%
%   ADet        The coupling coefficient of each detector.
%
%   Debug       OPTIONAL: Print out debugging info.
%
%   Calls: Hlm3DScatPntJacZB
%
%   Bugs: none known.

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
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: DPDWFDJacZB.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 17:45:29  rjg
%  Initial Revision for PMI 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PhiScat, PhiInc] = DPDWFDJacZB(CompVol, mu_sp, mu_a, delMu_a, v, f, ...
    pSrc, ASrc, pDet, ADet, Debug)

if nargin < 11
    Debug = 0;
end

%%
%%  Compute the PDE parameters from the media parameters
%%
D = v / (3 * (mu_sp + mu_a));
disp('Using David''s defn: D = v / (3 * mu_sp)')
D = v / (3 * (mu_sp));

ksqBkgnd = -v / D * mu_a  + j * 2*pi*f / D;
kBkgnd = sqrt(ksqBkgnd);

if Debug
    fprintf('D = %e\n', D);
    fprintf('Re{kBkgnd} = %f cm^-1\n', real(kBkgnd));
    fprintf('Im{kBkgnd} = %f cm^-1\n', imag(kBkgnd));
end

del_ksq = -v * delMu_a / D;
maxDelta = 1E-2;
nIterMax = 10000;

[PhiScat PhiInc nIter] = Hlm3DScatPntJacZB(CompVol, del_ksq, ksqBkgnd, ...
    pSrc, -v/D * ASrc, pDet, ADet, maxDelta, nIterMax, Debug);
