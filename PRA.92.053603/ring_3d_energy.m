%   Box size
L = 26;
Lz = 26;
%   Grid size
N=128;
Nz=128;
%   Initialize grid
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
k = [ (0:N/2)*2*pi/L -(N/2-1:-1:1)*2*pi/L];
kz = [ (0:Nz/2)*2*pi/Lz -(Nz/2-1:-1:1)*2*pi/Lz];
[KX,KY,KZ] = meshgrid(k,k,kz);
kk = ((KX.^2+KY.^2+KZ.^2)/2);

NE = (100);
g = 0.0276;
NN0=(1.e6); % number of particles
sr = 1;
s = 1;
omega = 1;

vm = 34;
vz = 18;
alpha = 0.8;
beta = 0.2;
m = 1;

re = [r(N/2+1:end)];
rze = r;
%%
[X,Z] = meshgrid(r,rz);
V = ((92-(0.5*Z.^2 * 9 + 0.5*X.^2 + 1.5*vm*beta^(2*m)*X.^2.^m.*exp(-m*(beta^2*X.^2 - 1)) + 1.5*vz*exp(-alpha^2*Z.^2)))).*tanh(sqrt((X).^2)/1).^2;
C=contour(X,Z,V,[3 3]);
%%
% plot(C(1,2:361),C(2,2:361));
% hold all;
% plot(C(1,364:401),C(2,364:401));
% plot(C(1,404:end),C(2,404:end));
z = C(2,2:222);
rr = C(1,2:222);
theta = 0:pi/20:pi*1.5;
xx = bsxfun(@times,rr',sin(theta));
yy = bsxfun(@times,rr',cos(theta));
zz = repmat(z',1,length(theta));
surf(xx,yy,zz,'FaceAlpha',0.6,'LineStyle','none','FaceLighting','gouraud',...
    'FaceColor',[0 1 0])
hold all
plot3(-rr, 0*rr, z, 'Color', [0 0.3 0], 'LineWidth', 1);
plot3(0*rr, rr, z, 'Color', [0 0.3 0], 'LineWidth', 1);
z = C(2,446:482);
rr = C(1,446:482);
theta = 0:pi/20:pi*1.5;
xx = bsxfun(@times,rr',sin(theta));
yy = bsxfun(@times,rr',cos(theta));
zz = repmat(z',1,length(theta));
surf(xx,yy,zz,'FaceAlpha',0.7,'LineStyle','none','FaceLighting','gouraud',...
    'FaceColor',[0 1 1])
light('Position',[-0.112527930568828 0.794552995758402 0.59667663082548]);
plot3([0 0], [0 0], [-5 5], 'Color', [0 0 1], 'LineWidth', 2);
theta = 0:pi/20:pi*2;
plot3(6*cos(theta), 6*sin(theta), 0*sin(theta), 'Color', [1 0 0], 'LineWidth', 2);
zlim([-10 10]);
view([27.5 18]);
plot3(-rr, 0*rr, z, 'Color', [0 0.3 0], 'LineWidth', 1);
plot3(0*rr, rr, z, 'Color', [0 0.3 0], 'LineWidth', 1);

phi = (0:300)/300*2*pi;
x = (6+1.7*cos(10*phi)).*cos(phi);
y = (6+1.7*cos(10*phi)).*sin(phi);
z = 0.8*sin(10*phi);
plot3(x, y, z, 'Color', [1 1 1], 'LineWidth', 1);

phi = (0:200)/200*2*pi;
x = (6+1.7*cos(10*phi)).*cos(phi);
y = (6+1.7*cos(10*phi)).*sin(phi);
z = 0.8*sin(10*phi);
u = gradient(x);
v = gradient(y);
w = gradient(z);
%quiver3(x,y,z,u,v,w,0);
arrow([x(1:10:200)' y(1:10:200)' z(1:10:200)'], [x(1:10:200)'+u(1:10:200)'/10 y(1:10:200)'+v(1:10:200)'/10 z(1:10:200)'+w(1:10:200)'/10], 4,...
    'FaceColor', [1 1 1], 'EdgeColor', [1 1 1], 'FaceLighting', 'none', 'LineWidth', 1);
%%
slice=zeros(N,Nz);
slice_ring=zeros(N,Nz);
[X,Y,Z] = meshgrid(r,r,rz);
RS = X.^2 + Y.^2;
V = 0.5*Z.^2 * 9 + 0.5*RS + vm*beta^(2*m)*RS.^m.*exp(-m*(beta^2*RS - 1)) + vz*exp(-alpha^2*Z.^2);
load('phi.mat'); %loading the initial state

xxi = zeros(length(re),Nz);


for ix=1:length(re)
    for iy=1:Nz
        xi = sqrt(0.5/g)/abs(phi(N/2,ix+N/2,iy));
        xxi(ix,iy)=xi;
    end
end
 xi=min(min(xxi));
 NN = sum(sum(sum(abs(phi.^2))));
 EEbase = sum(sum(sum(abs(conj(phi).*ifftn(kk.*fftn(phi)))  + (V+g*0.5*abs(phi.^2)).*abs(phi.^2) )))/NN;

%% line placement variations
EE=zeros(length(re),1);
EE2=zeros(length(re),1);
% r0 = 4.5;
% phi = phi.*tanh(sqrt((sqrt(RS)-r0).^2+(Z).^2)/xi).^sr.*exp(-1i*sr*atan2(sqrt(RS)-r0,Z));
for ix=1:length(re)
    phis = phi.*tanh(sqrt((X-re(ix)).^2+Y.^2)/xi).^s.*exp(-1i*s*atan2(X - re(ix),Y));
    phiv = abs(phis.^2);
    NN = sum(sum(sum(phiv)));
        EE(ix) = sum(sum(sum(abs(conj(phis).*ifftn(kk.*fftn(phis)))  + (V+g*0.5*phiv).*phiv )))/NN - EEbase;
        EE2(ix) =  real(sum(sum(sum( 1i*omega/(12*h)*conj(phis).*(X.*(-circshift(phis,2)+8*circshift(phis,1)-8*circshift(phis,-1)+circshift(phis,-2)) - ...
            Y.*(-circshift(phis,[0 2])+8*circshift(phis,[0 1])-8*circshift(phis,[0 -1])+circshift(phis,[0 -2]))) ))))/NN;
end
%% ring placement variations
EE=zeros(length(re),length(rze));
EE2=zeros(length(re),length(rze));
%  r0 = 0.2;
%  phi = phi.*tanh(sqrt((X-r0).^2+Y.^2)/xi).^s.*exp(-1i*s*atan2(X - r0,Y));
tic
for ix=1:length(re)
    for iz=1:length(rze)
%        phiv = abs(phi).^2.*tanh(sqrt(Y.^2+(X-re(ix)).^2)./xxi(ix)).^2;
%        phiv = abs(phi).^2.*tanh(sqrt((sqrt(RS)-re(ix)).^2+(Z-rze(iz)).^2)/xi).^2;
        phis = phi.*tanh(sqrt((sqrt(RS)-re(ix)).^2+(Z-rze(iz)).^2)/xi).^sr.*exp(-1i*sr*atan2(sqrt(RS)-re(ix),Z-rze(iz)));
        phiv = abs(phis.*conj(phis));
        NN = sum(sum(sum(phiv)));
%         EE(ix,iz) = 1/NN*0.5*sum(sum(sum( phiv./(Y.^2+(X-re(ix)).^2) )))*h*h*hz;
%         EE(ix,iz) = 1/NN*0.5*sum(sum(sum( phiv./((sqrt(RS)-re(ix)).^2+(Z-rze(iz)).^2) )))*h*h*hz;
        EE(ix,iz) = sum(sum(sum(abs(conj(phis).*ifftn(kk.*fftn(phis))) + (V+g*0.5*phiv).*phiv)))/NN - EEbase;
%        EE2(ix,iz) = real(sum(sum(sum( omega/(12*h)*conj(phis).*(X.*(-circshift(phis,2)+8*circshift(phis,1)-8*circshift(phis,-1)+circshift(phis,-2)) - ...
%            Y.*(-circshift(phis,[0 2])+8*circshift(phis,[0 1])-8*circshift(phis,[0 -1])+circshift(phis,[0 -2]))) ))))/NN;
%         phiv = abs(phi).^2.*tanh(sqrt(Y.^2+(X-re(ix)).^2)./min(xxi)).^2;
%        phiv = abs(phi).^2.*tanh(sqrt(Y.^2+(X-re(ix)).^2)./0.173785050707567).^2;
%        NN = sum(sum(sum(phiv)))*h*h*hz;
%        EE2(ix,iz) = 1/NN*0.5*sum(sum(sum( phiv./(Y.^2+(X-re(ix)).^2) )))*h*h*hz;
        %EE(N-ix+1,iz) = EE(ix,iz);
        %EE(ix,Nz-iz+1) = EE(ix,iz);
        %EE(N-ix+1,Nz-iz+1) = EE(ix,iz);
    end
end
toc
%plot(rze',EE(1,:)'/norm0);
% plot(re',EE(:,N/2)'/norm0);
%hold all;
%plot(re',xxi);
%% plotting density with energy contour
 figure;
 hold all;
 slice = squeeze(phi(:,N/2,:));
 imagesc(r,rz,abs(squeeze(phi(:,N/2,:))').^2);
 contour(re',rze',EE'/max(max(EE)),'Color',[1 1 1],'LevelStep',0.055);
%   contour(re',rze',EE2'/max(max(EE2)),'Color',[1 1 1],'LevelStep',0.05);
xlim([0 12]);
ylim([-12 12]);
colorbar;
%% saving cross-sections
% AAA = [r(N/2+1:end)' slice(N/2+1:end,N/2).^2];
% save density_z_0.dat  AAA -ASCII
% AAA = [r' slice(N/2+1,:)'.^2];
% save density_r_0.dat  AAA -ASCII
% AAA = [re' EE(:,NE+1)/10^6];
% save energy_z_0.dat  AAA -ASCII
% AAA = [rze' EE(1,:)'/10^6];
% save energy_r_0.dat  AAA -ASCII
