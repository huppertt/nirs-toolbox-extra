% function ds=medfiltPMI(ds,flag)
%
% medfiltPMI filtering ds.Recon.Mua with 3*3*3 median filtering
% method, and put the result in ds.Recon.MuaFilt
% usage:
%
% input:
%	ds: PMI data structure (structure to be presessed: ds.Recon.Mua
% output:
%	ds.Recon.MuaFilt

function ds=medfiltPMI(ds,flag)


if(nargin==2 & flag==1)

temp=reshape(ds.Recon.Musp(:,1),[size(ds.Inv.CompVol.Y,2),...
      size(ds.Inv.CompVol.X,2),size(ds.Inv.CompVol.Z,2)]);
temp=medfilt3(temp);
ds.Recon.MuspFilt(:,1)=reshape(temp,size(ds.Recon.Mua,1),1);
if(size(ds.Inv.Lambda,2)==2)
  temp=reshape(ds.Recon.Musp(:,2),[size(ds.Inv.CompVol.Y,2),...
      size(ds.Inv.CompVol.X,2),size(ds.Inv.CompVol.Z,2)]);
  temp=medfilt3(temp);
  ds.Recon.MuspFilt(:,2)=reshape(temp,size(ds.Recon.Mua,1),1);
end;

else

temp=reshape(ds.Recon.Mua(:,1),[size(ds.Inv.CompVol.Y,2),...
      size(ds.Inv.CompVol.X,2),size(ds.Inv.CompVol.Z,2)]);
temp=medfilt3(temp);
ds.Recon.MuaFilt(:,1)=reshape(temp,size(ds.Recon.Mua,1),1);
if(size(ds.Inv.Lambda,2)==2)
  temp=reshape(ds.Recon.Mua(:,2),[size(ds.Inv.CompVol.Y,2),...
      size(ds.Inv.CompVol.X,2),size(ds.Inv.CompVol.Z,2)]);
  temp=medfilt3(temp);
  ds.Recon.MuaFilt(:,2)=reshape(temp,size(ds.Recon.Mua,1),1);
end;

end;


