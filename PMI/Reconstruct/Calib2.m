function F=calib2(Mu,ds,k,l,Wavelength,idxFrame);
% Mu is the unknows we are solving for Mu_sp, Mu_a
%
%
%
%
%
%Mu_s_prime=Mu(1);
%Mu_a=Mu(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/05/25 13:14:47 $
%
%  $Revision: 1.1.1.1 $
%
%  $Log: Calib2.m,v $
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.3  2000/01/10 00:14:14  dboas
%  Storing the source and detector lists for use by other functions
%
%  Revision 1.2  1999/12/03 13:59:56  dboas
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

%
% g is the measurement data
%
g=reshape(ds.data.Raw(:,idxFrame,Wavelength),length(ds.data.Det.Pos),length(ds.data.Src.Pos))';
g_std=reshape(ds.data.Raw_std(:,idxFrame,Wavelength),length(ds.data.Det.Pos),length(ds.data.Src.Pos))';

if Mu(2) < 0
   F=4e50;
   return
end
   

%
% f is the calculated fluence at the detector using appro. model
%
ds.Fwd.Mu_sp(Wavelength) = Mu(1);
ds.Fwd.Mu_a(Wavelength) = Mu(2);
ds=genMeasData(ds);
f = ds.Fwd.PhiInc(:,Wavelength);

if length(f)==(length(ds.Fwd.Det.Pos)*length(ds.Fwd.Src.Pos))
   % MEASLIST IS UNIFORM
   f = reshape(f,length(ds.Fwd.Det.Pos),length(ds.Fwd.Src.Pos))';
   fd = f;
   gd = g;
   
   a=sum(g./(fd*diag((g./fd)'*(f(:,k)./gd(:,k)))),2);
	b=sum(g./(fd'*diag((g./fd)*(f(l,:)./gd(l,:))'))',1);
	L=a*b;
	inv_SD=(sum(sum(g.*f.*L))/sum(sum((L.*f).^2)));

	var = reshape(ds.data.Raw_std(:,idxFrame,Wavelength),length(ds.Fwd.Src.Pos),length(ds.Fwd.Det.Pos)).^2;
	sumand = (L.*f*inv_SD-g).^2;
	sumand = sumand;% ./ var;
	F=sum(sum(( sumand )));
else
   % MEASLIST IS NOT UNIFORM
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
   [c F] = calib_mex( gg, gg_std, ff );
end


%Pos_Source=squeeze(ds.data.Src.Pos(:,:,Wavelength));
%Pos_Detector=squeeze(ds.data.Det.Pos(:,:,Wavelength));
%f=PhiIncSR(Mu_a,Mu_s_prime,ds.data.speed/ds.data.idxRefr(idxFrame),0,ds.data.idxRefr(idxFrame),Pos_Source, Pos_Detector);
%f=reshape(f,16,9)';


disp([num2str(Mu) '  ' num2str(F)])