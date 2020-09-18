function dod = hmrMotionCorrectRLOESS(dod, t, span)
% Meryem Yucel, Oct, 2017
% span = 0.02 (default)
if span>0
    for i=1:size(dod,2)
    dod(:,i)=smooth(t,dod(:,i),span,'rloess');
    end
elseif span==-1
    return
end