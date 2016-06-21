h=0.15;
hz=0.2;
N=256;
Nz=32;
r=h*((1:N)-N/2-0.5);
rz=hz*((1:Nz)-Nz/2-0.5);
NE = 100;
g = 0.0276;
vm = 28;
vz = 10;
alpha = 0.8;
beta = 0.4;
m = 1;
%re = 0.1*(0:NE-1);
%rze = 0.1*(-NE:NE);
re = [r(N/2+1:end)];
rze = r;
%%
EE=zeros(length(re),1);
EE2=zeros(length(re),1);
slice=zeros(N,Nz);
slice_ring=zeros(N,Nz);
[X,Y,Z] = meshgrid(r,r,rz);
load('phi256.mat'); %loading the initial state

%% calculate energy
s = 1;
z0=2;
r0=5.5;
xi=0.12;
xxi = zeros(length(re),1);
EE=zeros(length(re),1);
NN0=6.0e5; % number of particles
tic
% norm=0;
% for ix2=N/2:N
%     for iz2=1:Nz
%         rr = sqrt(r(ix2).^2+r(N/2).^2);
%         norm = norm + abs(slice(ix2,iz2)).^2*rr;
%     end
% end
% norm0 = norm * 2*pi*h*hz;
for ix=1:length(re)

        xi = 5.155/abs(phi(N/2,ix+N/2,Nz/2));
        %xi = 4.255/abs(slice(N/2,N/2));
        xxi(ix)=xi;
end
for ix=1:length(re)

        phiv = abs(phi).^2.*tanh(sqrt(Y.^2+(X-re(ix)).^2)./xxi(ix)).^2;
        NN = sum(sum(sum(phiv)))*h*h*hz;
        EE(ix) = 1/NN*0.5*sum(sum(sum( phiv./(Y.^2+(X-re(ix)).^2) )))*h*h*hz;
        
%         phiv = abs(phi).^2.*tanh(sqrt(Y.^2+(X-re(ix)).^2)./min(xxi)).^2;
        phiv = abs(phi).^2.*tanh(sqrt(Y.^2+(X-re(ix)).^2)./0.173785050707567).^2;
        NN = sum(sum(sum(phiv)))*h*h*hz;
        EE2(ix) = 1/NN*0.5*sum(sum(sum( phiv./(Y.^2+(X-re(ix)).^2) )))*h*h*hz;
        %EE(N-ix+1,iz) = EE(ix,iz);
        %EE(ix,Nz-iz+1) = EE(ix,iz);
        %EE(N-ix+1,Nz-iz+1) = EE(ix,iz);

end
toc
%plot(rze',EE(1,:)'/norm0);
% plot(re',EE(:,N/2)'/norm0);
%hold all;
%plot(re',xxi);
%% plotting density with energy contour
% figure;
% hold all;
% image(r,rz,abs(slice').^2*64/38^2);
% contour(X,Z,EE'/max(max(EE))*max(max(slice))^2,'LineColor',[0 0 0],'LevelStep',80);
% xlim([0 10]);
% ylim([-10 10]);
% colorbar;
%% saving cross-sections
% AAA = [r(N/2+1:end)' slice(N/2+1:end,N/2).^2];
% save density_z_0.dat  AAA -ASCII
% AAA = [r' slice(N/2+1,:)'.^2];
% save density_r_0.dat  AAA -ASCII
% AAA = [re' EE(:,NE+1)/10^6];
% save energy_z_0.dat  AAA -ASCII
% AAA = [rze' EE(1,:)'/10^6];
% save energy_r_0.dat  AAA -ASCII
