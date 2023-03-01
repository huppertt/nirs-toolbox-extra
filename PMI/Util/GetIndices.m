function [index] = GetIndices(j,k,l,J,K,L)

% This function recieves the desired j,k,l indices and 
% returns the corresponding i indices in a jiffy.

js = length(j); 			
ks = length(k); 			
ls = length(l);

je = zeros(js,ks,ls);	
ke = zeros(js,ks,ls);	
le = zeros(js,ks,ls);

for a = 1:js
   je(a,:,:) = j(a);
end;

for a = 1:ks
   ke(:,a,:) = k(a);
end;

for a = 1:ls
   le(:,:,a) = l(a);
end;

ie = je+(ke-1)*J+(le-1)*J*K;

index = reshape(ie,js*ks*ls,1);

