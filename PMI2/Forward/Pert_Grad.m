% PERT_GRAD  Calcuate the gradient of a field.
%
% pert_grad = Pert_Grad(Pert, pmiModel, calc_mua )
%
%   This function calculates the gradiant of the field and is
%   used internally.  NOT DESIGNED FOR USE BY OTHERS.
%
%   Pert     - describes the field
%
%   Medium   - describes the Medium
%
%   calc_mua - flag to indicate if we also are interested in
%              mua perturbations
%
% Returns:
%      Pert_Grad - the gardient of Pert
%
% Known Problems or Limitations:
%   Assumes a voxel grid generated with meshgrid();
%   Assumes single wavelength, frequency domain only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2003, David Boas, Dana Brooks, Rick Gaudette, 
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

function Pert_Grad = Pert_Grad(Pert, Medium, calc_mua)

[X, Y, Z, dV] = sampleVolume(Medium.CompVol);

X = unique(X);    Nx = length(X);
Y = unique(Y);    Ny = length(Y);
Z = unique(Z);    Nz = length(Z);

Nsrcs = size(Pert,2);
Ndets = size(Pert,1) - Nx*Ny*Nz;

Pert_det = Pert(Ny*Nx*Nz + 1:size(Pert,1),:);
Pert_vox = reshape(Pert(1:Ny*Nx*Nz,:), Ny, Nx, Nz, Nsrcs);

Pert_x = zeros(size(Pert_vox));
Pert_y = zeros(size(Pert_vox));
Pert_z = zeros(size(Pert_vox));

if (Nx >= 2)
   if (Nx > 2)
      dX = X(3:Nx) - X(1:Nx-2);
      
      Pert_x(:,2:Nx-1,:,:) = ...
	  (Pert_vox(:,3:Nx,:,:) - Pert_vox(:,1:Nx-2,:,:)) ./ dX;
   end

   % End points
   Pert_x(:, 1,:,:) = ...
       (Pert_vox(:, 2,:,:) - Pert_vox(:,   1,:,:)) ./ (X(2) - X(1));
   Pert_x(:,Nx,:,:) = ...
       (Pert_vox(:,Nx,:,:) - Pert_vox(:,Nx-1,:,:)) ./ (X(Nx) - X(Nx-1));
end

if (Ny >= 2)
   if (Ny > 2)
      dY = Y(3:Ny) - Y(1:Ny-2);
      
      Pert_y(2:Ny-1,:,:,:) = ...
	  (Pert_vox(3:Ny,:,:,:) - Pert_vox(1:Ny-2,:,:,:)) ./ dY;
   end

   % End points
   Pert_y( 1,:,:,:) = ...
       (Pert_vox( 2,:,:,:) - Pert_vox(   1,:,:,:)) ./ (Y(2) - Y(1));
   Pert_y(Ny,:,:,:) = ...
       (Pert_vox(Ny,:,:,:) - Pert_vox(Ny-1,:,:,:)) ./ (Y(Ny) - Y(Ny-1));
end

if (Nz >= 2)
   if (Nz > 2)
      dZ = Z(3:Nz) - Z(1:Nz-2);
      
      Pert_z(:,:,2:Nz-1,:) = ...
	  (Pert_vox(:,:,3:Nz,:) - Pert_vox(:,:,1:Nz-2,:)) ./ dZ;
   end

   % End points
   Pert_z(:,:, 1,:) = ...
       (Pert_vox(:,:, 2,:) - Pert_vox(:,:,   1,:)) ./ (Z(2) - Z(1));
   Pert_z(:,:,Nz,:) = ...
       (Pert_vox(:,:,Nz,:) - Pert_vox(:,:,Nz-1,:)) ./ (Z(Nz) - Z(Nz-1));
end

Pert_vox = reshape(Pert_vox, Nx*Ny*Nz, Nsrcs);
Pert_det = reshape(Pert_det,    Ndets, Nsrcs);

Pert_x = reshape(Pert_x, Nx*Ny*Nz, Nsrcs);
Pert_y = reshape(Pert_y, Nx*Ny*Nz, Nsrcs);
Pert_z = reshape(Pert_z, Nx*Ny*Nz, Nsrcs);

if calc_mua
   Pert_Grad = [Pert_vox; Pert_x; Pert_y; Pert_z; Pert_det];
else
   Pert_Grad = [Pert_x; Pert_y; Pert_z; Pert_det];
end  

return;
