function C_sd=Calib_Coef(Mu,ds,Wavelength,idxFrame);
% Mu is the unknows we are solving for Mu_sp, Mu_a
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: Calib_Coef.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.3  2000/01/10 00:14:14  dboas
%  Storing the source and detector lists for use by other functions
%
%  Revision 1.2  1999/12/03 14:00:28  dboas
%  Modified to handle the case when measurements are not
%  made between every source and every detector.  This required
%  creation of a MEX file to do the nasty loops.
%
%  Revision 1.1  1999/11/18 14:19:34  tgaudett
%  Initial Calibration Routines
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mu_s_prime=Mu(1);
Mu_a=Mu(2);

%
% g is the measurment data
%
g=reshape(ds.data.Raw(:,idxFrame,Wavelength),length(ds.Fwd.Det.Pos),length(ds.Fwd.Src.Pos))';
g_std=reshape(ds.data.Raw_std(:,idxFrame,Wavelength),length(ds.Fwd.Det.Pos),length(ds.Fwd.Src.Pos))';

%
% f is the calculated fluence at the detector using appro. model
%
%Pos_Source=squeeze(ds.data.Src.Pos(:,:,Wavelength));
%Pos_Detector=squeeze(ds.data.Det.Pos(:,:,Wavelength));
%f=PhiIncSR(Mu_a,Mu_s_prime,ds.data.speed/ds.data.idxRefr(idxFrame),0,ds.data.idxRefr(idxFrame),Pos_Source, Pos_Detector);
f = ds.Fwd.PhiInc(:,Wavelength);

if length(f)==(length(ds.Fwd.Det.Pos)*length(ds.Fwd.Src.Pos))
	f=reshape(f,length(ds.Fwd.Det.Pos),length(ds.Fwd.Src.Pos))';

	ii=1;
	for l=1:length(ds.Fwd.Src.Pos)
   	for k=1:length(ds.Fwd.Det.Pos)
			a=sum(g./(f*diag((g./f)'*(f(:,k)./g(:,k)))),2);
			b=sum(g./(f'*diag((g./f)*(f(l,:)./g(l,:))'))',1);
			L=a*b;
			C_sd(ii) = sum(sum((L.*f).^2))/sum(sum(g.*f.*L));
      
	      ii=ii+1;
	   end
   end
   
else
   ff = zeros( length(ds.Fwd.Src.Pos), length(ds.Fwd.Det.Pos) );
   gg = zeros( length(ds.Fwd.Src.Pos), length(ds.Fwd.Det.Pos) );
   gg_std = zeros( length(ds.Fwd.Src.Pos), length(ds.Fwd.Det.Pos) );
   a = find(ds.Fwd.MeasList(:,4)==Wavelength);
   b = ds.Fwd.MeasList(a,:);
   for ii=1:length(a)
      ff( b(ii,1), b(ii,2) ) = f(ii);
      gg( b(ii,1), b(ii,2) ) = g( b(ii,1), b(ii,2) );
      gg_std( b(ii,1), b(ii,2) ) = g_std( b(ii,1), b(ii,2) );
   end
   [C_sd F] = calib_mex( gg, gg_std, ff );
   C_sd = reshape(C_sd,length(ds.Fwd.Src.Pos)*length(ds.Fwd.Det.Pos),1);
end


