% Pert_Grad(Pert, pmiModel, calc_mua )
%
%    This function calculates the gradiant of the field.  This
%    function is used by the routines HlmFullBorn_slab,
%    HlmFullBorn_ZB, and HlmFullBorn_NB for calculating the
%    perturbation matrix for the Nth Born approximation.
%
%    Pert - describes the field
%
%    pmiModel - describes the geometry of the field
%
%    calc_mua - flag to indicate if we also are interested in mua perturbations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%  $Author: dboas $
%
%  $Date: 2000/09/06 12:22:03 $
%
%  $Revision: 1.1 $
%
%  $Log: Pert_Grad.m,v $
%  Revision 1.1  2000/09/06 12:22:03  dboas
%  Initial revision.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Pert_Grad = Pert_Grad(Pert, pmiModel, calc_mua)

Nx = length(pmiModel.CompVol.X);
XStep = pmiModel.CompVol.XStep;
Ny = length(pmiModel.CompVol.Y);
YStep = pmiModel.CompVol.YStep;
Nz = length(pmiModel.CompVol.Z);
ZStep = pmiModel.CompVol.ZStep;

Nsrcs = size(Pert,2);
Ndets = size(Pert,1)-Nx*Ny*Nz;

Pert_det = Pert(Ny*Nx*Nz+1:size(Pert,1),:);
Pert_vox = reshape(Pert(1:Ny*Nx*Nz,:),[Ny, Nx, Nz, Nsrcs]);
Pert_x = zeros(size(Pert_vox));
Pert_y = zeros(size(Pert_vox));
Pert_z = zeros(size(Pert_vox));

if Nx>2
  Pert_x(:,2:Nx-1,:,:) = (Pert_vox(:,3:Nx,:,:) - Pert_vox(:,1: ...
				  Nx-2,:,:)) / (2*XStep);
end
if Nx>=2
  Pert_x(:,1,:,:) = (Pert_vox(:,2,:,:) - Pert_vox(:,1, ...
						  :,:)) / (XStep);
  Pert_x(:,Nx,:,:) = (Pert_vox(:,Nx,:,:) - Pert_vox(:,Nx-1, ...
						    :,:)) / (XStep);
end

if Ny>2
  Pert_y(2:Ny-1,:,:,:) = (Pert_vox(3:Ny,:,:,:) - Pert_vox(1:Ny-2,:, ...
				  :,:)) / (2*YStep);
end
if Ny>=2
  Pert_y(1,:,:,:) = (Pert_vox(2,:,:,:) - Pert_vox(1,:, ...
						  :,:)) / (YStep);
  Pert_y(Ny,:,:,:) = (Pert_vox(Ny,:,:,:) - Pert_vox(Ny-1,:, ...
						    :,:)) / (YStep);
end


if Nz>2
  Pert_z(:,:,2:Nz-1,:) = (Pert_vox(:,:,3:Nz,:) - Pert_vox(:,:, ...
				  1:Nz-2,:)) / (2*ZStep);
end
if Nz>=2
  Pert_z(:,:,1,:) = (Pert_vox(:,:,2,:) - Pert_vox(:,:, ...
						  1,:)) / (ZStep);
  Pert_z(:,:,Nz,:) = (Pert_vox(:,:,Nz,:) - Pert_vox(:,:, ...
						    Nz-1,:)) / (ZStep);
end


Pert_vox = reshape(Pert_vox,[Nx*Ny*Nz,Nsrcs]);
Pert_det = reshape(Pert_det,[Ndets,Nsrcs]);
Pert_x = reshape(Pert_x,[Nx*Ny*Nz,Nsrcs]);
Pert_y = reshape(Pert_y,[Nx*Ny*Nz,Nsrcs]);
Pert_z = reshape(Pert_z,[Nx*Ny*Nz,Nsrcs]);

if calc_mua
  Pert_Grad = [Pert_vox; Pert_x; Pert_y; Pert_z; Pert_det];
else
  Pert_Grad = [Pert_x; Pert_y; Pert_z; Pert_det];
end  