function [tIncCh, tInc] = hmrtInc_baselineshift_Ch(dod, t)

% Sahar Jahani, October 2017
iqr = 1.5;
fs = abs(1/(t(2)-t(1)));
for ww = 1:size(dod,2)
    clearvars -except ww tInc d ml fs t dod tMotion SD dodorg p iqr SNR_Thre SNR_Thre2 FrameSize_sec
    s1 =  dod(:,ww) ;
    
    %% Detecting motions
    fs=round(fs);
    [s1,ylpf] = hmrBandpassFilt( s1, fs, 0, 2);
    grad = conv(s1,[-1,0,1]);
    
    quants = quantile(grad,[.25 .50 .75]);  % compute quantiles
    IQR1 = quants(3)-quants(1);  % compute interquartile range
    High = quants(3)+IQR1*iqr;
    Low = quants(1)-IQR1*iqr;
    
    
    z=0;
    for i=1:length(dod)
        if ((grad(i)>High) || (grad(i)<Low))
            z=z+1; M(z)=i;
        end
    end

    extend = 12*fs;
    s11=repmat(s1(1,:),extend,1);s12=repmat(s1(end,:),extend,1);
    s1=[s11;s1;s12]; % extending signal for motion detection purpose (12 sec from each edge)
    
    fs=1/(t(2)-t(1));
    t=(0:(1/fs):(length(s1)/fs))';
    t=t(1:length(s1),1);
    
    %% Baseline shift motion detection
    
    if exist('M','var')
        M=M+extend;
        sig=ones(length(s1),1);
        for i=1:length(M)
            sig(M(i),:)=0;
        end
        
        %%% finding the location of the spikes or baseline shifts
        temp=(diff(sig));c1=0;c2=0;
        c=0;
        for i=1:length(s1)-1
            if temp(i)~=0
                c=c+1;
                pik(c)=i;
            end
        end
        temp2=diff(pik);
        
        c1=0;c2=0;
        for i=1:length(s1)-1
            if (temp(i)==1)
                c1=c1+1;
                meanpL(c1)=mean(s1(i),1);
            end
            
            if (temp(i)==-1)
                c2=c2+1;
                meanpH(c2)=mean(s1(i),1);
            end
        end
        
        motionkind=abs(meanpH-meanpL);
        
        %% finding the baseline shifts by comparing motion amplitudes with heart rate amplitude
        stemp=s1;
        [s1,ylpf] = hmrBandpassFilt( stemp, fs, 0, 2 );
        snoise2=stemp;
        zz=0;tt=1;
        for i=1:length(s1)-1
            if (sig(i)==1)
                zz=zz+1;
                sigtemp{1,tt}(1,zz)=s1(i);
                sigtempnoise{1,tt}(1,zz)=snoise2(i);
                if ((sig(i)==1) && (sig(i+1)==0))
                    tt=tt+1;
                    zz=0;
                end
            end
        end
        Nthre=round(0.5*fs);ssttdd=0;
        for i=1:tt
            tempo=sigtemp{1,i};
            if length(tempo)>Nthre
                for l=1:length(tempo)-Nthre
                    tempo2(l)=(abs(tempo(l+Nthre)-tempo(l)));
                end
            end
            ssttdd=[ssttdd tempo2];
            clear tempo2
            tempo2=[];
        end
        
        thrshld=quantile(ssttdd,[0.5]);
        pointS=(find(temp<0));
        pointE=(find(temp>0));
        countnoise=0;
        for ks=1:length(sigtempnoise)
            if (length(sigtempnoise{1,ks})>3*fs)
                countnoise=countnoise+1;
                dmean = mean(sigtempnoise{1,ks},2);
                dstd = std(sigtempnoise{1,ks},[],2);
                SNR_Thresh(countnoise,1)=abs(dmean)./dstd;
            end
        end
        
        SNR_Thre(1,ww)=mean(SNR_Thresh(2:end-1,1));
        
        
        sig2=ones(length(s1),1);
        
        for i=1:length(pointS)
            if motionkind(i)>thrshld
                sig2(pointS(i):pointE(i),:)=0;
            end
            % % % % % % % % % % % % % % % % % %
            % spline on long duration spikes  %
            % % % % % % % % % % % % % % % % % %
            
            if (((pointE(i)-pointS(i))> (0.1*fs))&&((pointE(i)-pointS(i))< (0.49999*fs)));
                sig2(pointS(i):pointE(i),:)=0;
            end
            if (pointE(i)-pointS(i))> (fs);
                sig2(pointS(i):pointE(i),:)=0;
            end
        end
        clear pointS
        clear pointE
        clear sig
        clear temp
        
        tInc(:,ww)=sig2;
    else
        tInc(:,ww)=ones(length(t),1);
    end
end

%% Extracting noisy channels from baseline-shift motion correction precedure
lent=size(dod,2)/2;
SNRvalue=3;
for ww=1:(size(dod,2)/2)-1
    if ((SNR_Thre(1,ww)<SNRvalue) && (SNR_Thre(1,ww+lent)<SNRvalue))
        tInc(:,ww+lent)=ones(length(t),1);
        tInc(:,ww)=ones(length(t),1);
    else if ((SNR_Thre(1,ww)>SNRvalue) && (SNR_Thre(1,ww+lent)<SNRvalue))
            tInc(:,ww+lent)=tInc(:,ww);
        else if ((SNR_Thre(1,ww)<SNRvalue) && (SNR_Thre(1,ww+lent)>SNRvalue))
                tInc(:,ww)=tInc(:,ww+lent);
            end
        end
    end
end

tIncCh=tInc(extend+1:end-extend,:);
tInc=tInc(extend+1:end-extend,:);
tIncall=tInc(:,1);
for kk=2:size(tInc,2)
  tIncall= tIncall.*tInc(:,kk);    
end
tInc=tIncall;
