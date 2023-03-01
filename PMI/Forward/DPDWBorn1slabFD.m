%DPDWBorn1slabFD    Forward DPDW model using 1st born for a slab medium,
%			but the 3pt fluence is calculated using Finite
%			difference.
%
%   [A Phi_Inc] = DPDWBorn1slabFD(Model, MeasList )
%
%   A           The forward matrix relating the contribution from each voxel to
%               a specific source-detector pair.  Each row is for a different
%               source detector pair.
%
%   Phi_Inc     The incident response at each detector from each source
%               in a column vector.  The same combination pattern as the
%               columns of A.
%
%   Model       The PMI Model structure contain the following fields: CompVol,
%               Mu_sp, Mu_a, v, idxRefr, and f.
%
%   MeasList    The measurement list for this simulation.  All frequencies and
%               wavelengths must be the same.
%
%   The true z boundary is assumed to be at z=0.
%
%   Calls: FDslab.m
%
%   Bugs: 
%        Assumes transmission measurements.
%

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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, Phi_Inc] = DPDWBorn1slabFD(Model, MeasList )

%%
%%  Get the frequency and wavelength index from MeasList
%%
idxFreq = MeasList(1,3);
if ~all(idxFreq == MeasList(:,3))
    error('all frequency indices must be the same in MeasList')
end
idxLambda = MeasList(1,4);
if ~all(idxLambda == MeasList(:,4))
    error('all wavelength indices must be the same in MeasList')
end


%%
%%  Compute the PDE parameters from the media parameters
%%
D = Model.v(idxLambda) / (3 * Model.Mu_sp(idxLambda));
k = sqrt(-Model.v(idxLambda) * Model.Mu_a(idxLambda) / D + ...
    j * 2 * pi * Model.ModFreq(idxFreq) * 1E6 / D);

%%
%%  Calculate the incident 2pt field
%%
[G, xfd, yfd, zfd] = FDslab( Model.Method.FD.xmin, ...
	    Model.Method.FD.xmax, ...
	    Model.Method.FD.ymin, ...
	    Model.Method.FD.ymax, ...
	    Model.Method.FD.zmin, ...
	    Model.Method.FD.zmax, ...
	    Model.Method.FD.step, ...
	    Model.Mu_sp(idxLambda), ...
	    Model.Mu_a(idxLambda) );

nx = size(G,1);
ny = size(G,2);
nz = size(G,3);
nxy = nx*ny;
Gd = G(:,:,nz:-1:1);
h = Model.Method.FD.step;

[Gx, Gy, Gz] = gradient( G, h );
[Gdx, Gdy, Gdz] = gradient( Gd, h );

%%
%%  Calculate the Born approximation and incident fields
%%
if strcmpi(Model.Method.Type,'Born') | strcmpi(Model.Method.Type,'Rytov')

  nMeas = size(MeasList,1);
  idxA = 1;

  if isfield(Model.Method,'ObjVec_mua')
    calc_mua = Model.Method.ObjVec_mua;
  else
    calc_mua = 0;
  end
  if isfield(Model.Method,'ObjVec_musp')
    calc_musp = Model.Method.ObjVec_musp;
  else
    calc_musp = 0;
  end

  for idxMeas = 1:nMeas
    
    xs = Model.Src.Pos( MeasList(idxMeas,1), 1);
    ys = Model.Src.Pos( MeasList(idxMeas,1), 2);
    xg = floor( (Model.CompVol.X - (Model.Method.FD.xmin+xs))/h );
    yg = floor( (Model.CompVol.Y - (Model.Method.FD.ymin+ys))/h ) + 1;
    zg = floor( (Model.CompVol.Z - (Model.Method.FD.zmin))/h );
    [xm, ym, zm] = meshgrid( xg, yg, zg );
    Slist = zm * nxy + xm * ny + ym;

    xd = Model.Det.Pos( MeasList(idxMeas,2), 1);
    yd = Model.Det.Pos( MeasList(idxMeas,2), 2);
    xg = floor( (Model.CompVol.X - (Model.Method.FD.xmin+xd))/h );
    yg = floor( (Model.CompVol.Y - (Model.Method.FD.ymin+yd))/h ) + 1;
    zg = floor( (Model.CompVol.Z - (Model.Method.FD.zmin))/h );
    [xm, ym, zm] = meshgrid( xg, yg, zg );
    Dlist = zm * nxy + xm * ny + ym;

    if calc_mua & ~calc_musp
      A(idxA,:) = G(Slist(:))' .* Gd(Dlist(:))';
    elseif ~calc_mua & calc_musp
      A(idxA,:) = Gx(Slist(:))' .* Gdx(Dlist(:))' + ...
                Gy(Slist(:))' .* Gdy(Dlist(:))' + ...
                Gz(Slist(:))' .* Gdz(Dlist(:))' ;
    elseif calc_mua & calc_musp
      A(idxA,:) = [G(Slist(:))' .* Gd(Dlist(:))' ...
                [Gx(Slist(:))' .* Gdx(Dlist(:))' + ...
                Gy(Slist(:))' .* Gdy(Dlist(:))' + ...
                Gz(Slist(:))' .* Gdz(Dlist(:))'] ] ;
    end

    xsd = floor( (xd-xs - Model.Method.FD.xmin) / h ) + 1;
    ysd = floor( (yd-ys - Model.Method.FD.ymin) / h ) + 1;
    Phi_Inc( idxA, 1 ) = G( (nz-1)*nxy + (xsd-1)*ny + ysd );

    idxA = idxA + 1;

  end

%  [A Phi_Inc] = Hlm3ptBorn1slab(Model, k, MeasList, zBnd, Debug);
  Phi_Inc = Model.v(idxLambda) / D * Phi_Inc;
  A = (Model.v(idxLambda) / D) * A;

elseif strcmpi(Model.Method.Type,'FullBorn') | ...
      strcmpi(Model.Method.Type,'BornN') | ...
      strcmpi(Model.Method.Type,'ExtBorn') | ...
      strcmpi(Model.Method.Type,'ExtBornN')

  error( 'Does not work for Method = slabFD' );

end

%%
%%  If Method = 'Rytov' then normalize the matrix by PhiInc
%%
if strcmp(Model.Method.Type,'Rytov')
  A = A ./ (Phi_Inc * ones(1,size(A,2)));
end
	


      
      