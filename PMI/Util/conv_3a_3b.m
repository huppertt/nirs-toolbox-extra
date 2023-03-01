%CONV_3A_3B     Convert a PMI structure from 3.0a to 3.0b
%
%   pmi = conv_3A_3B(pmi)
%
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
%  $Log: conv_3a_3b.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.2  1999/11/15 18:12:21  rjg
%  Added defaults for MeasList and Lambda
%
%  Revision 3.1  1999/11/10 19:01:08  rjg
%  Force all medium parameters to have the same dimension as the number of
%  wavelengths.
%  Fixed a bug in generating the number of Normals for the sources, it used to
%  generate the number of detectors.
%  Add MeasList and Lamba fields if they don't exist.  MeasList defaults to full
%  and Lamda defaults to 780 nm.
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pmi3b = conv_3a_3b(ds3a)

pmi3b = ds3a;
    
if isfield(ds3a, 'Fwd')
    if isfield(pmi3b.Fwd, 'nu')
        pmi3b.Fwd.v = pmi3b.Fwd.nu;
        pmi3b.Fwd = rmfield(pmi3b.Fwd, 'nu');
    end

    %%
    %%  Extend medium parameters to vectors if necessary for multiple wavelengths
    %%
    if isfield(pmi3b.Fwd, 'Lambda')
        nLambda = length(pmi3b.Fwd.Lambda);
    else
        nLambda = length(pmi3b.Fwd.Mu_a);
        pmi3b.Fwd.Lambda = 780;
        disp('INFO: Set Fwd.Lambda field to default value of 780 nm.');
    end
    if length(pmi3b.Fwd.idxRefr) ~= nLambda
        pmi3b.Fwd.idxRefr = repmat(pmi3b.Fwd.idxRefr, 1, nLambda);
    end
    if length(pmi3b.Fwd.g) ~= nLambda
        pmi3b.Fwd.g = repmat(pmi3b.Fwd.g, 1, nLambda);
    end
    if length(pmi3b.Fwd.v) ~= nLambda
        pmi3b.Fwd.v = repmat(pmi3b.Fwd.v, 1, nLambda);
    end
    if length(pmi3b.Fwd.Mu_s) ~= nLambda
        pmi3b.Fwd.Mu_s = repmat(pmi3b.Fwd.Mu_s, 1, nLambda);
    end
    if length(pmi3b.Fwd.Mu_sp) ~= nLambda
        pmi3b.Fwd.Mu_sp = repmat(pmi3b.Fwd.Mu_sp, 1, nLambda);
    end

    pmi3b.Fwd.Src = ds3a.Fwd.SrcPos;
    pmi3b.Fwd.Det = ds3a.Fwd.DetPos;
    pmi3b.Fwd = rmfield(pmi3b.Fwd, 'SrcPos');
    pmi3b.Fwd = rmfield(pmi3b.Fwd, 'DetPos');

    %%
    %%  Handle the old source detector plane idea
    %%
    if strcmp(ds3a.Fwd.SrcPos.Type, 'uniform')
        nSrc = length(ds3a.Fwd.SrcPos.X) * length(ds3a.Fwd.SrcPos.Y) * ...
            length(ds3a.Fwd.SrcPos.Z);
    else
        nSrc = length(ds3a.Fwd.SrcPos.Pos);
    end
    if strcmp(ds3a.Fwd.DetPos.Type, 'uniform')
        nDet = length(ds3a.Fwd.DetPos.X) * length(ds3a.Fwd.DetPos.Y) * ...
            length(ds3a.Fwd.DetPos.Z);
    else
        nDet = length(ds3a.Fwd.DetPos.Pos);
    end

    if length(ds3a.Fwd.SrcAmp) > 1
        nXYSrc = nSrc / length(ds3a.Fwd.SrcPos.Z);
        amp = 10 ^ (ds3a.Fwd.SrcAmp ./ 20)
        amp = repmat(amp(:)', nXYSrc, 1);
        pmi3b.Fwd.Src.Amplitude = amp(:);
    else
        pmi3b.Fwd.Src.Amplitude = repmat(10 ^ (ds3a.Fwd.SrcAmp ./ 20), nSrc, 1);
    end
    %%
    %%  Assume diffuse medium is -z domain
    %%
    pmi3b.Fwd.Src.Normal = repmat([0 0 -1], nSrc, 1);
    pmi3b.Fwd.Det.Normal = repmat([0 0 -1], nDet, 1);
    pmi3b.Fwd.Det.Amplitude = ones(nDet, 1);
    pmi3b.Fwd.Src.NA = ones(nSrc, 1);
    pmi3b.Fwd.Det.NA = ones(nDet);
    pmi3b.Fwd = rmfield(pmi3b.Fwd, 'SrcAmp');
    pmi3b.Fwd = rmfield(pmi3b.Fwd, 'SensorError');

    %%
    %%  Expand the medium parameters if necessary and create the wavelength
    %%  parameter
    %%
    
    %%
    %%  Boundary description
    %%
    switch lower(ds3a.Fwd.Boundary)
     case {'extrapolated', 'semi', 'semi-infinite'}
      pmi3b.Fwd = rmfield(pmi3b.Fwd, 'Boundary');
      pmi3b.Fwd.Boundary.Geometry = 'Semi-infinite';
     case {'infinite', 'inf'}
      pmi3b.Fwd = rmfield(pmi3b.Fwd, 'Boundary');
      pmi3b.Fwd.Boundary.Geometry = 'Infinite';
     
     otherwise
      error(['Unknown boundary condition: ' ds3a.Fwd.Boundary.Geometry]);
    end

    %%
    %%  Add the MeasList table to the data structure if it does not exist.
    %%
    if ~isfield(pmi3b.Fwd, 'MeasList')
        pmi3b.Fwd.MeasList = FullMeasList(nSrc, nDet, ...
            length(pmi3b.Fwd.ModFreq), length(pmi3b.Fwd.Lambda));
    end
    
end

if isfield(ds3a, 'Inv')
    if isfield(pmi3b.Inv, 'nu')
        pmi3b.Inv.v = pmi3b.Inv.nu;
        pmi3b.Inv = rmfield(pmi3b.Inv, 'nu');
    end
    %%
    %%  Extend medium parameters to vectors if necessary for multiple wavelengths
    %%
    if isfield(pmi3b.Inv, 'Lambda')
        nLambda = length(pmi3b.Inv.Lambda);
    else
        nLambda = length(pmi3b.Inv.Mu_a);
        pmi3b.Inv.Lambda = 780;
        disp('INFO: Set Inv.Lambda field to default value of 780 nm.');
    end
    if length(pmi3b.Inv.idxRefr) ~= nLambda
        pmi3b.Inv.idxRefr = repmat(pmi3b.Inv.idxRefr, 1, nLambda);
    end
    if length(pmi3b.Inv.g) ~= nLambda
        pmi3b.Inv.g = repmat(pmi3b.Inv.g, 1, nLambda);
    end
    if length(pmi3b.Inv.v) ~= nLambda
        pmi3b.Inv.v = repmat(pmi3b.Inv.v, 1, nLambda);
    end
    if length(pmi3b.Inv.Mu_s) ~= nLambda
        pmi3b.Inv.Mu_s = repmat(pmi3b.Inv.Mu_s, 1, nLambda);
    end
    if length(pmi3b.Inv.Mu_sp) ~= nLambda
        pmi3b.Inv.Mu_sp = repmat(pmi3b.Inv.Mu_sp, 1, nLambda);
    end

    
    pmi3b.Inv.Src = ds3a.Inv.SrcPos;
    pmi3b.Inv.Det = ds3a.Inv.DetPos;
    pmi3b.Inv = rmfield(pmi3b.Inv, 'SrcPos');
    pmi3b.Inv = rmfield(pmi3b.Inv, 'DetPos');

    %%
    %%  Handle the old source detector plane idea
    %%
    if strcmp(ds3a.Inv.SrcPos.Type, 'uniform')
        nSrc = length(ds3a.Inv.SrcPos.X) * length(ds3a.Inv.SrcPos.Y) * ...
            length(ds3a.Inv.SrcPos.Z);
    else
        nSrc = length(ds3a.Inv.SrcPos.Pos);
    end
    if strcmp(ds3a.Inv.DetPos.Type, 'uniform')
        nDet = length(ds3a.Inv.DetPos.X) * length(ds3a.Inv.DetPos.Y) * ...
            length(ds3a.Inv.DetPos.Z);
    else
        nDet = length(ds3a.Inv.DetPos.Pos);
    end

    if length(ds3a.Inv.SrcAmp) > 1
        nXYSrc = nSrc / length(ds3a.Inv.SrcPos.Z);
        amp = 10 ^ (ds3a.Inv.SrcAmp ./ 20)
        amp = repmat(amp(:)', nXYSrc, 1);
        pmi3b.Inv.Src.Amplitude = amp(:);
    else
        pmi3b.Inv.Src.Amplitude = repmat(10^(ds3a.Inv.SrcAmp./20),nSrc,1);
    end

    %%
    %%  Assume diffuse medium is -z domain
    %%
    pmi3b.Inv.Src.Normal = repmat([0 0 -1], nSrc, 1);
    pmi3b.Inv.Det.Normal = repmat([0 0 -1], nDet, 1);
    pmi3b.Inv.Det.Amplitude = ones(nDet, 1);
    pmi3b.Inv.Src.NA = ones(nSrc, 1);
    pmi3b.Inv.Det.NA = ones(nDet);
    pmi3b.Inv = rmfield(pmi3b.Inv, 'SrcAmp');
    pmi3b.Inv = rmfield(pmi3b.Inv, 'SensorError');

    %%
    %%  Boundary description
    %%
    switch lower(ds3a.Inv.Boundary)
     case {'extrapolated', 'semi', 'semi-infinite'}
      pmi3b.Inv = rmfield(pmi3b.Inv, 'Boundary');
      pmi3b.Inv.Boundary.Geometry = 'Semi-infinite';
     case {'infinite', 'inf'}
      pmi3b.Inv = rmfield(pmi3b.Inv, 'Boundary');
      pmi3b.Inv.Boundary.Geometry = 'Infinite';
     
     otherwise
      error(['Unknown boundary condition: ' ds3a.Inv.Boundary.Geometry]);
    end

    if isfield(ds3a.Inv, 'Method')
        pmi3b.Inv = rmfield(pmi3b.Inv, 'Method');
        pmi3b.Inv.Method.Type = ds3a.Inv.Method;
        if any(strcmp(pmi3b.Inv.Method.Type, 'Born'))
            pmi3b.Inv.Method.Order = 1;
        end
    end

    %%
    %%  Add the MeasList table to the data structure if it does not exist.
    %%
    if ~isfield(pmi3b.Inv, 'MeasList')
        pmi3b.Inv.MeasList = FullMeasList(nSrc, nDet, ...
            length(pmi3b.Inv.ModFreq), length(pmi3b.Inv.Lambda));
    end

end

if isfield(ds3a, 'Object')
    pmi3b = rmfield(pmi3b, 'Object');
    nSpheres = 0;
    if isfield(ds3a.Object, 'SphereCtr')
        nSpheres = size(ds3a.Object.SphereCtr, 1);
        for i = 1:nSpheres
            pmi3b.Object{i}.Type = 'sphere';
            pmi3b.Object{i}.Pos = ds3a.Object.SphereCtr(i,:);
            pmi3b.Object{i}.Radius = ds3a.Object.SphereRad(i);
            pmi3b.Object{i}.Mu_a = ds3a.Object.SphereDelta + ds3a.Fwd.Mu_a;
        end
    end
    if isfield(ds3a.Object, 'BlockCtr')
        nBlocks = size(ds3a.Object.BlockCtr, 1);
        for i = 1:nBlocks
            pmi3b.Object{i+nSpheres}.Type = 'block';
            pmi3b.Object{i+nSpheres}.Pos = ds3a.Object.BlockCtr(i,:);
            pmi3b.Object{i+nSpheres}.Radius = ds3a.Object.BlockDims(i,:);
            pmi3b.Object{i+nSpheres}.Mu_a = ds3a.Object.BlockDelta + ...
                ds3a.Fwd.Mu_a;
        end
    end
end
