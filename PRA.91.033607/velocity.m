% load('params.mat')
% load('slice_00322.mat')
h = 0.3;
N = 128;
r=h*((1:N)-N/2-0.5);
[X,Y] = meshgrid(r,r);
[sx,sy] = gradient(slice);
[sxc,syc] = gradient(conj(slice));
vx = -imag(slice.*sxc - conj(slice).*sx)./(abs(slice).^2).*(abs(slice)>1e-4)./(2*h);
vy = -imag(slice.*syc - conj(slice).*sy)./(abs(slice).^2).*(abs(slice)>1e-4)./(2*h);
[XX,YY] = meshgrid(r,r);
aXY = (atan2(YY,XX));
vn = vx.*cos(aXY) + vy.*sin(aXY);
vt = -vx.*sin(aXY) + vy.*cos(aXY);
%%

npoints = 600;

ang = linspace(180,360,npoints)/180*pi;
dphi = ang(2)-ang(1);
% 322
%  R = 14.49;
% R = 6.066;

% 305
% R = 14.76;
% R = 5.903;
R = 10.45;
xarr = R*cos(ang);
yarr = R*sin(ang);

distr = interp2(r,r,slice,xarr,yarr);
VN = interp2(r,r,vn,xarr,yarr); %/(sqrt(MU(322)/2));
VT = interp2(r,r,vt,xarr,yarr); %/(sqrt(MU(322)/2));

% sum(VT)*R*dphi/(2*pi)
% plot(ang(1:end-1)/pi*180,diff(interp2(r,r,angle(slice),xarr,yarr))/dphi/R);
% plot(ang/pi*180,interp2(r,r,angle(slice),xarr,yarr));


%%
startr = [ -6:-0.7:-14 6:0.7:14];
cphase = 2.5*2*pi*80*1.29366e-3*4;
starty = startr*sin(cphase);
startx = startr*cos(cphase);

imagesc(r,r,abs(slice.^2));
set(gca,'YDir','normal');
hold on;
s1 = streamline(X,Y,vx,vy,startx,starty);
s2 = streamline(X,Y,-vx,-vy,startx,starty);
set(s1,'Color',[1,1,1]);
set(s2,'Color',[1,1,1]);
set(s1,'LineWidth',1.5);
set(s2,'LineWidth',1.5);

%%

ang = linspace(0,360,npoints)/180*pi;
R = 14.76;
xarr = R*cos(ang);
yarr = R*sin(ang);
plot(xarr,yarr,'Color',[0.9,0.9,0.9],'LineWidth',1.5);
R = 10.45;
xarr = R*cos(ang);
yarr = R*sin(ang);
plot(xarr,yarr,'Color',[0.9,0.9,0.9],'LineWidth',1.5);
R = 5.903;
xarr = R*cos(ang);
yarr = R*sin(ang);
plot(xarr,yarr,'Color',[0.9,0.9,0.9],'LineWidth',1.5);