%STRIPRES       Strip the results from a PMI data structure.
%
%   pmi = stripres(pmi)
%
%   pmi         The PMI data structure.
%
%
%   STRIPRES will strip the results fields from the PMI data strucutre specfied.
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
%  $Log: stripres.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 3.0  1999/06/17 21:18:09  rjg
%  Initial PMI 3.0 revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ds = stripres(ds)

if isfield(ds, 'PhiTotal')
    ds = rmfield(ds, 'PhiTotal');
end

if isfield(ds, 'PhiTotalN')
    ds = rmfield(ds, 'PhiTotalN');
end

if isfield(ds, 'PhiTotalNw')
    ds = rmfield(ds, 'PhiTotalNw');
end
if isfield(ds, 'IPR')
    ds = rmfield(ds, 'IPR');
end

if isfield(ds, 'Fwd')
    if isfield(ds.Fwd, 'A');
        ds.Fwd = rmfield(ds.Fwd, 'A');
    end
    if isfield(ds.Fwd, 'PhiInc')
        ds.Fwd = rmfield(ds.Fwd , 'PhiInc');
    end
    if isfield(ds.Fwd, 'PhiScat')
        ds.Fwd = rmfield(ds.Fwd , 'PhiScat');
    end
    if isfield(ds.Fwd, 'PhiScatw')
        ds.Fwd = rmfield(ds.Fwd , 'PhiScatw');
    end
end

if isfield(ds, 'Inv')
    if isfield(ds.Inv, 'A');
        ds.Inv = rmfield(ds.Inv, 'A');
    end
    if isfield(ds.Inv, 'Aw');
        ds.Inv = rmfield(ds.Inv, 'Aw');
    end
    if isfield(ds.Inv, 'PhiInc')
        ds.Inv = rmfield(ds.Inv , 'PhiInc');
    end
    if isfield(ds.Inv, 'PhiScat')
        ds.Inv = rmfield(ds.Inv , 'PhiScat');
    end
    if isfield(ds.Inv, 'PhiScatw')
        ds.Inv = rmfield(ds.Inv , 'PhiScatw');
    end
end


if isfield(ds, 'Noise')
    if isfield(ds.Noise, 'w');
        ds.Noise = rmfield(ds.Noise, 'w');
    end
    if isfield(ds.Noise, 'SrcNoise');
        ds.Noise = rmfield(ds.Noise, 'SrcNoise');
    end
    if isfield(ds.Noise, 'DetNoise');
        ds.Noise = rmfield(ds.Noise, 'DetNoise');
    end
    if isfield(ds.Noise, 'ScatNoise');
        ds.Noise = rmfield(ds.Noise, 'ScatNoise');
    end
end
