%CONV_2_3       Convert a DPDW 2.0 structure to a PMI 3.0a strucutre
%
%   pmi = conv_2_3(dpdw)
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
%  $Log: conv_2_3.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ds3 = conv_2_3(ds2)

ds3.Version = 3.0;
ds3.Debug = 1;
ds3.Fwd.CompVol = ds2.CompVol;
ds3.Fwd.SrcPos = ds2.SrcPos;
ds3.Fwd.DetPos = ds2.DetPos;
ds3.Fwd.SrcAmp = ds2.SrcAmp;
ds3.Fwd.SensorError = ds2.SensorError;
ds3.Fwd.ModFreq = ds2.ModFreq;
ds3.Fwd.idxRefr = ds2.idxRefr;
ds3.Fwd.g = ds2.g;
ds3.Fwd.Mu_s = ds2.Mu_s;
ds3.Fwd.Mu_a = ds2.Mu_a;
ds3.Fwd.Mu_sp = ds2.Mu_sp;
ds3.Fwd.nu = ds2.nu;
ds3.Fwd.Boundary = ds2.Boundary;
switch ds2.DataSource
    case 'Born-1'
        ds3.Fwd.Method.Type = 'Born';
        ds3.Fwd.Method.Order = 1;

    case 'Rytov-1'
        ds3.Fwd.Method.Type = 'Rytov';
        ds3.Fwd.Method.Order = 1;
end

ds3.Inv.CompVol = ds2.CompVol;
ds3.Inv.SrcPos = ds2.SrcPos;
ds3.Inv.DetPos = ds2.DetPos;
ds3.Inv.SrcAmp = ds2.SrcAmp;
ds3.Inv.SensorError = ds2.SensorError;
ds3.Inv.ModFreq = ds2.ModFreq;
ds3.Inv.idxRefr = ds2.idxRefr;
ds3.Inv.g = ds2.g;
ds3.Inv.Mu_s = ds2.Mu_s;
ds3.Inv.Mu_a = ds2.Mu_a;
ds3.Inv.Mu_sp = ds2.Mu_sp;
ds3.Inv.nu = ds2.nu;
ds3.Inv.Boundary = ds2.Boundary;
switch ds2.DataSource
    case 'Born-1'
        ds3.Inv.Method = 'Born';
        ds3.Inv.BornOrder = 1;

    case 'Rytov-1'
        ds3.Inv.Method = 'Rytov';
        ds3.Inv.RytovOrder = 1;
end

ds3.Noise.SrcSNRflag = ds2.SrcSNRflag;
if ds2.SrcSNRflag
    ds3.Noise.SrcSNR = ds2.SrcSNR;
end
ds3.Noise.DetSNRflag = ds2.DetSNRflag;
if ds2.DetSNRflag
    ds3.Noise.DetSNR = ds2.DetSNR;
end
ds3.Noise.ScatSNRflag = ds2.ScatSNRflag;
if ds2.ScatSNRflag
    ds3.Noise.ScatSNR = ds2.ScatSNR;
end

if isfield(ds2, 'VecNormSNRflag')
    ds3.Noise.VecNormSNRflag = ds2.VecNormSNRflag;
    if ds2.VecNormSNRflag
        ds3.Noise.VecNormSNR = ds2.VecNormSNR;
    end
else
    ds3.Noise.VecNormSNRflag = 0;
end


if isfield(ds2, 'SphereCtr')
    ds3.Object.SphereCtr = ds2.SphereCtr;
    ds3.Object.SphereRad = ds2.SphereRad;
    ds3.Object.SphereDelta = ds2.SphereDelta;
end

if isfield(ds2, 'BlockCtr')
    ds3.Object.BlockCtr = ds2.BlockCtr;
    ds3.Object.BlockRad = ds2.BlockRad;
    ds3.Object.BlockDelta = ds2.BlockDelta;
end

ds3.Recon.ReconAlg = ds2.ReconAlg;
ds3.Visualize.VisPlane = ds2.VisPlane;
ds3.Visualize.Type = ds2.Visualize;
ds3.Visualize.PlaneIndices = ds2.PlaneIndices;
ds3.Visualize.LayoutVector = ds2.LayoutVector;
ds3.Visualize.CMap = ds2.CMap;
ds3.Visualize.CRange = ds2.CRange;
if isfield(ds2, 'nCLines')
    ds3.Visualize.nCLines = ds2.nCLines;
end
