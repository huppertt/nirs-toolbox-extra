%SHOWIMAGE        
%
%	fig = showimage(ds);
%   Calls: 
%
%   Bugs: none known.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2001/02/16 13:38:08 $
%
%  $Revision: 1.1 $
%  Initial revision
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fig = showimage(ds);

if isfield(ds.Recon,'Mua');
   Mu_a = reshape(ds.Recon.Mua,[size(ds.Inv.CompVol.Y,2),...
	    size(ds.Inv.CompVol.X,2), size(ds.Inv.CompVol.Z,2),...
	    size(ds.Inv.Lambda,2)]);
   maxMua= max(Mu_a(:));
   minMua= min(Mu_a(:));
   for c = 1:size(ds.Inv.Lambda,2)
       figure;
       nZ = size(ds.Inv.CompVol.Z,2);
       ngx=1; while ceil(nZ/ngx)>ngx
	 ngx=ngx+1;
       end
       if nZ<2;
	  imagesc(ds.Inv.CompVol.X,ds.Inv.CompVol.Y,...
		   squeeze(Mu_a(:,:,1,c)),[minMua maxMua]);
	  colorbar
	  title(...
	  sprintf('mu_a lambda=%d, Z=%5.3f',ds.Inv.Lambda(c),...
		ds.Inv.CompVol.Z(1)));
	else
           for cc = 1:nZ
	      subplot(ngx,ngx,cc)
	      imagesc(ds.Inv.CompVol.X,ds.Inv.CompVol.Y,...
		   squeeze(Mu_a(:,:,cc,c)),[minMua maxMua]);
	       colorbar
	       title(...
	       sprintf('mu_a lambda=%d, Z=%5.3f',ds.Inv.Lambda(c),...
		   ds.Inv.CompVol.Z(cc)));
	   end;
         end;
    end;
end;
if isfield(ds.Recon,'Musp');
   Mu_sp = reshape(ds.Recon.Musp,[size(ds.Inv.CompVol.Y,2),...
	    size(ds.Inv.CompVol.X,2), size(ds.Inv.CompVol.Z,2),...
	    size(ds.Inv.Lambda,2)]);
   maxMusp= max(Mu_sp(:));
   minMusp= min(Mu_sp(:));
   for c = 1:size(ds.Inv.Lambda,2)
       figure;
       nZ = size(ds.Inv.CompVol.Z,2);
       ngx=1; while ceil(nZ/ngx)>ngx
	 ngx=ngx+1;
       end
       if nZ<2;
	  imagesc(ds.Inv.CompVol.X,ds.Inv.CompVol.Y,...
		   squeeze(Mu_sp(:,:,1,c)),[minMusp maxMusp]);
	  colorbar
	  title(...
	  sprintf('mu_sp lambda=%d, Z=%5.3f',ds.Inv.Lambda(c),...
		ds.Inv.CompVol.Z(1)));
	else
           for cc = 1:nZ
	      subplot(ngx,ngx,cc)
	      imagesc(ds.Inv.CompVol.X,ds.Inv.CompVol.Y,...
		   squeeze(Mu_sp(:,:,cc,c)),[minMusp maxMusp]);
	       colorbar
	       title(...
	       sprintf('mu_sp lambda=%d, Z=%5.3f',ds.Inv.Lambda(c),...
		   ds.Inv.CompVol.Z(cc)));
	   end;
         end;
    end;
end;
