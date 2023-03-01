% function[SD,Phi,data,hdr,mask] = packset(basefile, threshold);

function[h] = showmask(mask)

if (isempty(mask))
   error('Mask cannot be empty list');
end

if (~iscell(mask))
   error('Not a cell array');
end

img = full(mask{1});

for k = 2:length(mask)
   ml = find(mask{k});
   img(ml) = img(ml) + 1;
end

h = figure;
imagesc(img);

return;

