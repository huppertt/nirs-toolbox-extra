% function[Phi,Measlist] = repack_image(data, dark, mask)

function[Phi,MeasList] = repack_image(data, dark, mask)

[n1,n2,nsrc,ntim] = size(data);
[n3,n4,ndet]      = size(mask);
[n5,n6]           = size(dark);

if ((n1 ~= n3) | (n2 ~= n4))
   error('data and mask images have different sizes');
end

if ((n1 ~= n5) | (n2 ~= n6))
   error('data and dark images have different sizes');
else
   dark = double(dark);
end

% Allocate memory

MeasList = zeros(ndet * nsrc * ntim, 9);
Phi      = zeros(ndet * nsrc * ntim, 1);

% Normalize mask to 1

msk = cell(1,ndet);

for idet = 1:ndet;
   msk{idet} = sparse(mask(:,:,idet) / sum(sum(mask(:,:,idet))));
end

% Turn image into measurements

imeas = 1;

for itim = 1:ntim;
   for isrc = 1:nsrc;
      img = double(data(:,:,isrc,itim)) - dark;

      for idet = 1:ndet;
         MeasList(imeas,1:6) = [ isrc idet 0 1 1 itim ];
%         Phi(imeas) = sum(sum( img .* mask(:,:,idet) ));
         Phi(imeas) = sum(sum( img .* msk{idet} ));
         imeas = imeas + 1;
      end;
   end;

   disp([ 'itim = ', num2str(itim) ]);
end

return;
