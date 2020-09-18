function dod = hmrMotionCorrectSG(dod,t, FrameSize_sec)
% Savitzky_Golay filter
% K: polynomial order (default = 3)
% FrameSize_sec: Frame size in sec (default = 10 sec)
% Meryem Yucel, Oct, 2017

K = 3;
fs = abs(round(1/(t(1)-t(2))));
FrameSize_sec = (FrameSize_sec * fs); 
if mod(FrameSize_sec,2)==0 
    FrameSize_sec = FrameSize_sec  + 1;
end

dod=sgolayfilt(dod,K,FrameSize_sec);