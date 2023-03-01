%geomDisplay       Display Geometry of the problem for talks
%
%   fig = geomDisplay(pmi)
%
%   pmi         The Photon Migration Imaging data structure.
%
%
%   Returns the figure number of the figure generated.
%
%   geomDisplay displays the geometry of the problem labeling sources,
%     detectors and objects in the correct orientation and size.  This also
%     helps in visualizing the problem.
%
%
%
%   Calls: none.
%
%   Bugs: It does not handle a Sphere all that well but it gives you and idea
%	 We also need to handle a MeasList for showing Sources and Detectors, 
%	this is not implimented yet.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  $Author: tgaudett $
%
%  $Date: 2000/10/17 15:16:27 $
%
%  $Revision: 1.3 $
%  Initial Revision for PMI 3.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function geomDisplay(pmi);
%
%  We need to first check to see if both the forward and inverse
%   models are setup and whether they are the same or not.
%	
Fwd_Flag 	= 0;
Inv_Flag 	= 0;
Object_Flag 	= 0;

if(isfield(pmi,'Fwd'))
	%  We have a Fwd Model
	Fwd = pmi.Fwd;
	Fwd_Flag = 1;
end;
if(isfield(pmi,'Fwd'))
	%  We have a Inv Model
	Inv = pmi.Inv;
	Inv_Flag = 1;
end;
if(isfield(pmi,'Object'))
	%  We have an Object
	Object = pmi.Object;
	Object_Flag = 1;
end;

%
%	Now we need to find out how big space needs to be.
%
minmax = getminmax(pmi);
figure;
if(Fwd_Flag)
	fig  = dispModel(Fwd,minmax);
	fig1 = dispObject(Object,Fwd,fig);
	title('Forward Model');
end;
figure;
if(Inv_Flag)
	fig2 = dispModel(Inv,minmax);
	fig3 = dispObject(Object,Inv,fig2);
	title('Inverse Model');
end;
 
fig = [fig,fig1,fig2,fig3];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Internal Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig = dispModel(Model,minmax);
%
%
%

axis(minmax);
hold on;
plot3(Model.Src.Pos(:,1),Model.Src.Pos(:,2),Model.Src.Pos(:,3),'ro',...
	Model.Det.Pos(:,1),Model.Det.Pos(:,2),Model.Det.Pos(:,3),'bo');
legend('Sources','Detectors');
%
% Draw mesh on Bottom.  If I fill the hole mesh it getts to busy.
% X-Y Plane Z = 0;
if Model.CompVol.Z(1) < 0
	zp = length(Model.CompVol.Z);
else
	zp = 1;
end;
for y=1:length(Model.CompVol.Y)-1
for x=1:length(Model.CompVol.X)-1
line([Model.CompVol.X(x),Model.CompVol.X(x+1)],...
	[Model.CompVol.Y(y),Model.CompVol.Y(y)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(x+1),Model.CompVol.X(x+1)],...
	[Model.CompVol.Y(y),Model.CompVol.Y(y+1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(x+1),Model.CompVol.X(x)],...
	[Model.CompVol.Y(y+1),Model.CompVol.Y(y+1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(x),Model.CompVol.X(x)],...
	[Model.CompVol.Y(y+1),Model.CompVol.Y(y)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');

end;
end;

%
%  Now just draw box for full CompVol
%
if zp==1
	zp = length(Model.CompVol.Z);
else
   zp = 1;
end;
	
xp = length(Model.CompVol.X);
yp = length(Model.CompVol.Y);

%
%  X-Y Plane Z = zp
%
line([Model.CompVol.X(1),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(1)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');

zp = length(Model.CompVol.Z);

%
%  X-Z Plane Y = 0
%
line([Model.CompVol.X(x),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(1)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(1)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(1)],'Color','b');


%
%  X-Z Plane Y = yp
%
line([Model.CompVol.X(x),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(1)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(1)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(1)],'Color','b');

%
%  Y-Z Plane X = 0
%
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(1)],'Color','b');
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(1),Model.CompVol.X(1)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(1)],'Color','b');


%
%  Y-Z Plane X = xp
%
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(1)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(yp)],...
	[Model.CompVol.Z(1),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(yp),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(zp)],'Color','b');
line([Model.CompVol.X(xp),Model.CompVol.X(xp)],...
	[Model.CompVol.Y(1),Model.CompVol.Y(1)],...
	[Model.CompVol.Z(zp),Model.CompVol.Z(1)],'Color','b');





fig = gcf;
hold off;
return;
function fig=dispObject(Object,Model,fig);
%
%
%

figure(fig);

%
%	Display Object
%
color = ['blue ';'red  ';'green'];
for c=1:length(Object)
	if strcmp(Object{c}.Type,'Block')
		Vp=GenBlock(Model.CompVol,Object{c}.Pos,Object{c}.Dims,1);
	else
		Vp=GenSphere(Model.CompVol,Object{c}.Pos,Object{c}.Radius,1);
	end;
	[X,Y,Z]=meshgrid(Model.CompVol.X,Model.CompVol.Y,...
			Model.CompVol.Z);
	idxblock = find(Vp==1);
	for nblock = 1:length(idxblock)
		patchblock(X(idxblock(nblock)),Y(idxblock(nblock)),...
			Z(idxblock(nblock)),Model.CompVol.XStep,...
			Model.CompVol.YStep,...
			Model.CompVol.ZStep,color(c,:));
	end;		
	camlight; %camlight(-80,-10);
	lighting phong
end;


return;

function minmax = getminmax(pmi);
Fwd_Flag 	= 0;
Inv_Flag 	= 0;
Object_Flag 	= 0;

if(isfield(pmi,'Fwd'))
	%  We have a Fwd Model
	Fwd = pmi.Fwd;
	Fwd_Flag = 1;
end;
if(isfield(pmi,'Fwd'))
	%  We have a Inv Model
	Inv = pmi.Inv;
	Inv_Flag = 1;
end;
if(isfield(pmi,'Object'))
	%  We have an Object
	Object = pmi.Object;
	Object_Flag = 1;
end;

max_x = 0;
max_y = 0;
max_z = 0;
min_x = 0;
min_y = 0;
min_z = 0;
if(Fwd_Flag)
	max_x = max([Fwd.Src.Pos(:,1);Fwd.Det.Pos(:,1);Fwd.CompVol.X';max_x]);
	max_y = max([Fwd.Src.Pos(:,2);Fwd.Det.Pos(:,2);Fwd.CompVol.Y';max_y]);
	max_z = max([Fwd.Src.Pos(:,3);Fwd.Det.Pos(:,3);Fwd.CompVol.Z';max_z]);
	min_x = min([Fwd.Src.Pos(:,1);Fwd.Det.Pos(:,1);Fwd.CompVol.X';min_x]);
	min_y = min([Fwd.Src.Pos(:,2);Fwd.Det.Pos(:,2);Fwd.CompVol.Y';min_y]);
	min_z = min([Fwd.Src.Pos(:,3);Fwd.Det.Pos(:,3);Fwd.CompVol.Z';min_z]);
end;
if(Inv_Flag)
	max_x = max([Inv.Src.Pos(:,1);Inv.Det.Pos(:,1);Inv.CompVol.X';max_x]);
	max_y = max([Inv.Src.Pos(:,2);Inv.Det.Pos(:,2);Inv.CompVol.Y';max_y]);
	max_z = max([Inv.Src.Pos(:,3);Inv.Det.Pos(:,3);Inv.CompVol.Z';max_z]);
	min_x = min([Inv.Src.Pos(:,1);Inv.Det.Pos(:,1);Inv.CompVol.X';min_x]);
	min_y = min([Inv.Src.Pos(:,2);Inv.Det.Pos(:,2);Inv.CompVol.Y';min_y]);
	min_z = min([Inv.Src.Pos(:,3);Inv.Det.Pos(:,3);Inv.CompVol.Z';min_z]);
end;
if(Object_Flag)
	for c = 1:length(Object)
		if strcmp(Object{c}.Type,'Block')
		    max_x = max([Object{c}.Pos(:,1)+Object{c}.Dims(1)/2;...
				     max_x]);
		    max_y = max([Object{c}.Pos(:,2)+Object{c}.Dims(2)/2;...
				     max_y]);
		    max_z = max([Object{c}.Pos(:,3)+Object{c}.Dims(3)/2;...
				     max_z]);
		    min_x = min([Object{c}.Pos(:,1)-Object{c}.Dims(1)/2;...
				     min_x]);
		    min_y = min([Object{c}.Pos(:,2)-Object{c}.Dims(2)/2;...
				     min_y]);
		    min_z = min([Object{c}.Pos(:,3)-Object{c}.Dims(3)/2;...
				     min_z]);
		else
		    max_x = max([Object{c}.Pos(:,1)+Object{c}.Radius(1)/2;...
				     max_x]);
		    max_y = max([Object{c}.Pos(:,2)+Object{c}.Radius(1)/2;...
				     max_y]);
		    max_z = max([Object{c}.Pos(:,3)+Object{c}.Radius(1)/2;...
				     max_z]);
		
		    min_x = min([Object{c}.Pos(:,1)-Object{c}.Radius(1)/2;...
				     min_x]);
		    min_y = min([Object{c}.Pos(:,2)-Object{c}.Radius(1)/2;...
				     min_y]);
		    min_z = min([Object{c}.Pos(:,3)-Object{c}.Radius(1)/2;...
				     min_z]);
		end;
	end;

end;
minmax = [min_x,max_x,min_y,max_y,min_z,max_z];
return

function p=patchblock(x,y,z,xstep,ystep,zstep,color)

xstep = xstep/2;
ystep = ystep/2;
zstep = zstep/2;

%
%	Side 1
%
xp = [x-xstep,x+xstep,x+xstep,x-xstep,x-xstep];
yp = [y-ystep,y-ystep,y+ystep,y+ystep,y-ystep];
zp = [z-zstep,z-zstep,z-zstep,z-zstep,z-zstep];
p(1)=patch(xp,yp,zp,color);

%
%	Side 2
%
zp = [z+zstep,z+zstep,z+zstep,z+zstep,z+zstep];
p(2)=patch(xp,yp,zp,color);

%
%	Side 3
%
xp = [x-xstep,x-xstep,x-xstep,x-xstep,x-xstep];
yp = [y-ystep,y-ystep,y+ystep,y+ystep,y-ystep];
zp = [z-zstep,z+zstep,z+zstep,z-zstep,z-zstep];
p(3)=patch(xp,yp,zp,color);

%
%	Side 4
%
xp = [x+xstep,x+xstep,x+xstep,x+xstep,x+xstep];
p(4)=patch(xp,yp,zp,color);

%
%	Side 5
%
xp = [x-xstep,x-xstep,x+xstep,x+xstep,x-xstep];
yp = [y-ystep,y-ystep,y-ystep,y-ystep,y-ystep];
zp = [z-zstep,z+zstep,z+zstep,z-zstep,z-zstep];
p(5)=patch(xp,yp,zp,color);

%
%	Side 6
%
yp = [y+ystep,y+ystep,y+ystep,y+ystep,y+ystep];
p(6)=patch(xp,yp,zp,color);


