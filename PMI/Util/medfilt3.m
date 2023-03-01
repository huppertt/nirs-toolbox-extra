function tomoImg=medfilt3(data)
%medfilt3 is a 3 dimensional median filtering program, filter pattern is 3*3*3
%input: data: 3 dimentional X*Y*Z matrix
%output: tomoImg: median filtered data
if(size(data,1)<3|size(data,2)<3|size(data,3)<3)
   error('the data matrix is too small. Each dimension must be larger than 3')
end;
sizeX=size(data,1);
sizeY=size(data,2);
sizeZ=size(data,3);
tempImg=zeros(sizeX+2,sizeY+2,sizeZ+2);
tempImg(2:sizeX+1,2:sizeY+1,2:sizeZ+1)=data;
cntrX=1;
cntrY=1;
cntrZ=1;
%filtPat=data(cntrX-1:cntrX+1,cntrY-1:cntrY+1,cntrZ-1:cntrZ+1);
for(ZZ=0:sizeZ-1)
   for(YY=0:sizeY-1)
      for(XX=0:sizeX-1)
         tomoImg(cntrX+XX,cntrY+YY,cntrZ+ZZ)= median(reshape...
            (tempImg(cntrX+XX:cntrX+XX+2,cntrY+YY:cntrY+YY+2,cntrZ+ZZ:cntrZ+ZZ+2),1,27));
      end;
   end;
end;


