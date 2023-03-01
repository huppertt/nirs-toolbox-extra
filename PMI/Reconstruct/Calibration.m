%function ds=Calibration(ds,idxCalibrateFrame,param,debug);
%
%	
% User must supply an intial guess for every frame.  This is
% initial guess is provided in ds.data.Mu_sp and ds.data.Mu_a.
%
%
%
% Internal Functions
%		RUM
%			PhiIncSR  Incident Fluence at the Recievers form sources using Semi-Infinite Medium. Extrap Boundary
%		RCM
%		FCM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: dboas $
%
%  $Date: 2000/07/27 15:08:46 $
%
%  $Revision: 1.3 $
%
%  $Log: Calibration.m,v $
%  Revision 1.3  2000/07/27 15:08:46  dboas
%  Many fixes to handle the new data-structure arrangement for raw data and MODEL
%  data.
%
%  Revision 1.2  2000/06/27 14:36:17  dboas
%  Modifactions needed to handle special situations when small subsets of the
%  MeasList are used for the calibration.
%
%  Revision 1.1.1.1  2000/05/25 13:14:47  dboas
%  initial
%
%  Revision 1.2  1999/12/03 14:04:32  dboas
%  Modified to handle calibration of the raw_std data as well as the raw data.
%
%  Revision 1.1  1999/11/18 14:19:34  tgaudett
%  Initial Calibration Routines
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ds=Calibration(ds,idxCalibrateFrame,param,debug);

  if ~exist('param')
    param = [];
  end
  if ~isfield(param,'lambda')
    param.lambda = 1e-1;
  end
  if ~isfield(param,'lambdaInc')
    param.lambdaInc = 2.0;
  end
  if ~isfield(param,'lambdaDec')
    param.lambdaDec = 1.5;
  end
  if ~isfield(param,'lambdaMax')
    param.lambdaMax = 1e4;
  end
  if ~isfield(param,'verbose')
    param.verbose = 0;
  end


if(strcmp(ds.data.CalibrationMethod,'RefUncorMeas'));
%  ds.data.PhiIncSR = ds.Inv.P.PhiInc;
  % Loop over all Frames they want to work on
  oldCalFrame=0;
  for nWorkingOnFrame=1:length(idxCalibrateFrame)
    iFrame=idxCalibrateFrame(nWorkingOnFrame);
    nCalFrame=find(ds.data.CalFlag(1:iFrame)==1);
    nCalFrame=nCalFrame(length(nCalFrame));
    % Calculate Corrected Data
    [ds]=CalibBase(ds,iFrame,nCalFrame,oldCalFrame==nCalFrame);
    oldCalFrame = nCalFrame;   
  end;
  
elseif(strcmp(ds.data.CalibrationMethod,'RefUncorMeasOpt'));
  % Loop over all Frames they want to work on
  oldCalFrame=0;
  for nWorkingOnFrame=1:length(idxCalibrateFrame)
    iFrame=idxCalibrateFrame(nWorkingOnFrame);
    nCalFrame=find(ds.data.CalFlag(1:iFrame)==1);
    nCalFrame=nCalFrame(length(nCalFrame));
    % Calculate Corrected Data
    [ds]=CalibOpt(ds,iFrame,nCalFrame,oldCalFrame==nCalFrame,0,1,param,debug);
    oldCalFrame = nCalFrame;   
  end;
%  ds.data.PhiIncSR = ds.Inv.PhiInc;
  
elseif(strcmp(ds.data.CalibrationMethod,'FrmUncorMeasOpt'));
  % Loop over all Frames they want to work on
  for nWorkingOnFrame=1:length(idxCalibrateFrame)
    iFrame=idxCalibrateFrame(nWorkingOnFrame);
%    nCalFrame=find(ds.data.CalFlag(1:iFrame)==1);
%    nCalFrame=nCalFrame(length(nCalFrame));
    % Calculate Corrected Data
    [ds]=CalibOpt(ds,iFrame,iFrame,0,0,1,param,debug);
  end;
  %	   ds.data.PhiIncSR = ds.Fwd.PhiInc;

elseif(strcmp(ds.data.CalibrationMethod,'RefCorMeasOpt'));
  % Loop over all Frames they want to work on
  oldCalFrame=0;
  for nWorkingOnFrame=1:length(idxCalibrateFrame)
    iFrame=idxCalibrateFrame(nWorkingOnFrame);
    nCalFrame=find(ds.data.CalFlag(1:iFrame)==1);
    nCalFrame=nCalFrame(length(nCalFrame));
    % Calculate Corrected Data
    [ds]=CalibOpt(ds,iFrame,nCalFrame,oldCalFrame==nCalFrame,1,1,param,debug);
    oldCalFrame = nCalFrame;   
  end;
%  ds.data.PhiIncSR = ds.Fwd.PhiInc;
  
elseif(strcmp(ds.data.CalibrationMethod,'FrmCorMeasOpt'));
  % Loop over all Frames they want to work on
  for nWorkingOnFrame=1:length(idxCalibrateFrame)
    iFrame=idxCalibrateFrame(nWorkingOnFrame);
%    nCalFrame=find(ds.data.CalFlag(1:iFrame)==1);
%    nCalFrame=nCalFrame(length(nCalFrame));
    % Calculate Corrected Data
    [ds]=CalibOpt(ds,iFrame,iFrame,0,1,1,param,debug);
  end;

elseif(strcmp(ds.data.CalibrationMethod,'FrmCorMeas'));
  % Loop over all Frames they want to work on
  for nWorkingOnFrame=1:length(idxCalibrateFrame)
    iFrame=idxCalibrateFrame(nWorkingOnFrame);
    [ds]=CalibOpt(ds,iFrame,iFrame,0,1,0,param,debug);
  end;

else
  disp('WARNING: Calibration Method is unknown');

end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Internal Functions 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calibrate given a baseline measurement (nCalFrame)
% and the background optical properties.
% Do not consider the measurements to be correlated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ds]=CalibBase(ds,iFrame,nCalFrame,CalibrationFlag);

Nw=ds.data.nWavelengths;

%
% Calibrate only against the SDpairs specified in the Fwd Model.
%
if CalibrationFlag==0

  for idxLambda=1:Nw
    meas = ConvertModel2Complex(ds.Inv, ds.Inv.P(idxLambda).PhiInc);
    [raw raw_std] = RawDataSubset( ds.Inv, ds.data, idxLambda, nCalFrame );
    M(:,idxLambda) = raw ./ meas;
    badNumbers=find(M(:,idxLambda)==0);
    M(badNumbers,idxLambda)=1;
  end      
  ds.data.Coupling(:,nCalFrame)=M(:);
end

ds.data.Corrected(:,iFrame) = ds.data.Raw(:,iFrame) ./ ...
    ds.data.Coupling(:,nCalFrame);
ds.data.Corrected_std(:,iFrame) = ds.data.Raw_std(:,iFrame) ./ ...
    ds.data.Coupling(:,nCalFrame);



%for idxLambda=1:Nw   
%  [raw raw_std] = RawDataSubset( ds.Inv, ds.data, idxLambda, iFrame );
%  ds.data.Corrected(:,iFrame,idxLambda)= raw(:) ./ ds.data.Coupling(:,nCalFrame,idxLambda);
%  ds.data.Corrected_std(:,iFrame,idxLambda)= raw_std(:) ./ ds.data.Coupling(:,nCalFrame,idxLambda);
%end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calibrate a data set (nCalFrame) and fit
% for the background optical properties of that
% data set.  The determined optical properties
% are stored as are the coupling coefficients.
% These coupling coefficients are then used to correct
% the iFrames data sets.
%
% The calibration can be done either considering
% the correlation of the measurements or not.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ds]=CalibOpt(ds,iFrames,nCalFrame,CalibrationFlag, ...
		       CorrelationFlag, OptFlag, param,debug);
  
nWavelengths = length(ds.Inv.Lambda);

if( CalibrationFlag==0)
  [raw raw_std] = RawDataSubset( ds.Inv, ds.data, nCalFrame );

  clear temp_PhiTotalN;
  if isfield(ds.Inv,'PhiTotalN')
    temp_PhiTotalN = ds.Inv.PhiTotalN;
  end
  clear ds.Inv.PhiTotalN;
  
  if ds.Inv.ModFreq>0
    ds.Inv.PhiTotalN = [real(raw);imag(raw)];
  else
    ds.Inv.PhiTotalN = real(raw);
  end
  
  
  mu(1,:) = ds.data.Mu_sp(nCalFrame,:);
  mu(2,:) = ds.data.Mu_a(nCalFrame,:);
  if CorrelationFlag
    if OptFlag
      [fmu lsq CCs CCd PhiInc] = FitBackgroundCal( mu, ds, param );
    else
      fmu = mu;
      for idxLambda=1:length(ds.Inv.Lambda)
	foo = find(idxLambda==ds.Inv.MeasList(:,4));
	if ~isempty(foo)
	  [delta_mu lsq CCs(:,idxLambda) CCd(:,idxLambda) PhiInc(:,idxLambda)] = ...
	      FitBackgroundCal1( mu(:,idxLambda), ds.Inv, idxLambda, ...
				 param.lambda);
	  if param.verbose
	    lsq
	  end
	end
      end
    end
  else
    [fmu lsq CCs CCd PhiInc] = FitBackground( mu, ds );
  end

  clear ds.Inv.PhiTotalN;
  if exist('temp_PhiTotalN')
    ds.Inv.PhiTotalN = temp_PhiTotalN;
  end
  
  if debug > 0
    foos = sprintf( 'Optical Properties for Frame %d, [%5.2f\t%5.4f', nCalFrame, fmu(1,1), fmu(2,1) );
    for idxLambda = 2:nWavelengths
      foos = sprintf( '%s\t%5.2f\t%5.4f', foos, fmu(1,idxLambda), fmu(2,idxLambda) );
    end
    foos = sprintf( '%s]', foos );
    disp( foos );
  end
  
  ds.data.Mu_s(nCalFrame,:) = fmu(1,:);
  ds.data.Mu_sp(nCalFrame,:) = fmu(1,:);
  ds.data.Mu_a(nCalFrame,:) = fmu(2,:);

  ds.data.Coupling(:,nCalFrame) = ones(size(ds.data.Raw,1),1);
  ds.data.PhiInc(:,nCalFrame) = zeros(size(ds.data.Raw,1),1);
  ds.data.PhiTotalN(:,nCalFrame) = zeros(size(ds.data.Raw,1),1);

  for idx=1:size(ds.Inv.MeasList,1)
    idxLambda = ds.Inv.MeasList(idx,4);
    idxSrc = ds.Inv.MeasList(idx,1);
    idxDet = ds.Inv.MeasList(idx,2);

%    idxLambdaList = find( idxLambda==ds.data.MeasList(:,4) );
    idxList = find( idxSrc==ds.data.MeasList(:,1) & ...
		    idxDet==ds.data.MeasList(:,2) & ...
		    idxLambda==ds.data.MeasList(:,4) );
    if ~isempty(idxList)
      ds.data.Coupling(idxList, nCalFrame) = ...
	  CCs(idxSrc,idxLambda) * CCd(idxDet,idxLambda);
      
      idxFoo = find( idxLambda==ds.Inv.MeasList(:,4) );
      idxMeas = find( idxSrc==ds.Inv.MeasList(idxFoo,1) & ...
		      idxDet==ds.Inv.MeasList(idxFoo,2));
      ds.data.PhiInc(idxList,nCalFrame) = PhiInc(idxMeas, ...
						 idxLambda);
      ds.data.PhiTotalN(idxList,nCalFrame) = ds.data.Raw(idxList, ...
			nCalFrame) ./ ds.data.Coupling(idxList,nCalFrame);
    end
  end
end

for idxLambda = 1:nWavelengths
  ds.data.Mu_s(iFrames,idxLambda) = ds.data.Mu_s(nCalFrame,idxLambda);
  ds.data.Mu_sp(iFrames,idxLambda) = ds.data.Mu_sp(nCalFrame,idxLambda);
  ds.data.Mu_a(iFrames,idxLambda) = ds.data.Mu_a(nCalFrame,idxLambda);
%  ds.data.Coupling(iFrames,idxLambda) = ds.data.Coupling(nCalFrame,idxLambda);

  ds.data.Coupling(:,iFrames) = ds.data.Coupling(:,nCalFrame);
  ds.data.PhiInc(:,iFrames) = ds.data.PhiInc(:,nCalFrame);
  ds.data.PhiTotalN(:,iFrames) = ds.data.PhiTotalN(:,nCalFrame);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Determine the subset of Raw data specified by the Fwd Model
% MeasList
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raw, raw_std] = RawDataSubset( Model, RawData, idxLambda, ...
					 nCalFrame)



%
% Calibrate only against the SDpairs specified in the Fwd Model.
%
%nLambda = length(Model.Lambda);

%for idxLambda=1:nLambda
  a = find(Model.MeasList(:,4)==idxLambda);
%  if ~isempty(a)
    aa = Model.MeasList(a,:);
%    aa = Model.MeasList(:,:);
    bb = RawData.MeasList(:,:);
    index1 = zeros(length(aa),1);
    for ii=1:length(aa)
      index1(ii) = find( bb(:,1)==aa(ii,1) & ...
			 bb(:,2)==aa(ii,2) & ...
			 bb(:,3)==aa(ii,3) & ...
			 bb(:,4)==aa(ii,4) );
    end

    raw = squeeze(RawData.Raw(index1(:),nCalFrame));
    raw_std = squeeze(RawData.Raw_std(index1(:), ...
					 nCalFrame));
%  end
%end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert the Model data to complex
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meas = ConvertModel2Complex( Model, ModelData )

  if Model.ModFreq > 0
    nMeas = size(ModelData,1)/2;
    meas = ModelData(1:nMeas) + sqrt(-1) * ModelData(nMeas+1:2* ...
						     nMeas);
  else
    meas = ModelData;
  end
  
   
