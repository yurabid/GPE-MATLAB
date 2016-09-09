global tcoef Ub omegaBar initphase ub1;
%   Box size
xmax = 15;
zmax = 7.5;
%   Grid size
N=128;
Nz=32;
grid = grid3d(xmax,N,xmax,N,zmax,Nz);

hbar=1.054571800e-34;
mass=1.41922608e-25;
omegascale = 570*2*pi;
tcoef = 1/(omegascale); %time scale
rscale = sqrt(hbar/mass/omegascale);

Ub=1.60931; %*1.5
ub1 = 0.0;
omegaBar = 0.0;
initphase = 0;%-0.125*pi;

r0 = 8.85787;
omz = (0.526316);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + (sqrt(X.^2+Y.^2)-r0).^2);

task = GPEtask(grid,Vfun);
task.g = 0.160567;%/sqrt(2*pi/omz);
task.Ntotal = 5000; 
task.Vtd = @TDPot;


[phi, mu1] = task.groundstate_itp(1e-2,1e-6);
mu=grid.inner(phi,task.applyham(phi))/task.Ntotal; %mu1(end);
cs=sqrt(task.g*max(max(max(abs(task.init_state).^2)))/2)/r0;
%% config.figure
% figure;
% hold all;
% VV = task.getVtotal(0);
% plot(grid.x*rscale*1e6,VV(64,:,16)*570/1000);
% plot(grid.x*rscale*1e6,ub1*grid.x*570/1000);
% 
% xvals = grid.x(find(VV(64,:,16)<real(mu(end)+0.06*r0) & grid.x<0))*rscale*1e6;
% yvals = VV(64,find(VV(64,:,16)<real(mu(end)+0.06*r0) & grid.x<0),16)*570/1000;
% baseval = real(mu(end)+0.06*r0)*570/1000;
% aaa=patch([xvals flip(xvals)], [yvals yvals*0+baseval], 'red');
% % aaa.BaseLine.LineStyle = 'none';
% xvals = grid.x(find(VV(64,:,16)<real(mu(end)-0.06*r0) & grid.x>0))*rscale*1e6;
% yvals = VV(64,find(VV(64,:,16)<real(mu(end)-0.06*r0) & grid.x>0),16)*570/1000;
% baseval = real(mu(end)-0.06*r0)*570/1000;
% bbb=patch([xvals flip(xvals)], [yvals yvals*0+baseval], 'red');
% % bbb=area(grid.x(find(VV(64,:,16)<real(mu(end)-0.06*r0) & grid.x>0))*rscale*1e6,VV(64,find(VV(64,:,16)<real(mu(end)-0.06*r0) & grid.x>0),16)*570/1000,real(mu(end)-0.06*r0)*570/1000);
% % bbb.BaseLine.LineStyle = 'none';
% 
% plot(grid.x*rscale*1e6,0*grid.x + real(mu(end)-0.06*r0)*570/1000);
% plot(grid.x*rscale*1e6,0*grid.x + real(mu(end)+0.06*r0)*570/1000);

%%
nlev = 1;
spectr = zeros(60,nlev);

rr=grid.x(N/2+1:end);
Nr = length(rr);
[RR,ZZ] = meshgrid(rr,grid.z);
angles = atan2(grid.mesh.x,grid.mesh.y);
asangles = (1:256)/256*2*pi;
[ras,pas] = meshgrid(rr,asangles);
% xas = ras.*sin(pas);
% yas = ras.*cos(pas);
[xas,yas] = pol2cart(pas,ras);
Vrz = reshape(Vfun(RR,0,ZZ),Nr*Nz,1);
phirz = reshape(squeeze(phi(N/2,N/2+1:end,:)).',Nr*Nz,1);
lapr = findiff1_arb(rr,2,0)+spdiags(1./rr(:),0,Nr,Nr)*findiff1_arb(rr,1,0);
lapz = findiff1_arb(grid.z,2,0);
nlin = task.g*spdiags(abs(phirz).^2,0,Nr*Nz,Nr*Nz);


plot(rr,abs(phi(N/2,N/2+1:end,Nz/2).^2)./task.Ntotal);
% plot((0:40),(0:40)*cs*570/1000,'--');
opts.v0 = [phirz; phirz];
for i=7:7 %[1,2,10,20,30]
m=i-1;
figure;
% hold all;
h1 = -0.5*(kron(lapr,speye(Nz)) + kron(speye(Nr),lapz) - m^2*kron(spdiags(1./rr(:).^2,0,Nr,Nr),speye(Nz))) + spdiags(Vrz,0,Nr*Nz,Nr*Nz) + 2*nlin - (mu-0.00705)*speye(Nr*Nz);
mtot = [h1 nlin; -nlin -h1];
[vv,dd]=eigs(mtot,nlev,'sm',opts);
ddd=diag(dd);
spectr(i,:) = ddd(:);
% plot((1:nlev)*0+m,abs(real(ddd))*570/1000,'LineStyle','none','Marker','o','Color',[0 0 1],'MarkerFaceColor',[0 0 1],'MarkerSize',4);
% plot((1:nlev)*0+m,abs(imag(ddd)),'LineStyle','none','Marker','.','Color',[1 0 0]);
vdens = (abs(reshape(conj(phirz).*(vv(1:end/2,nlev)+conj(vv(end/2+1:end,nlev))),Nz,Nr)));
vres = reshape(vv(1:Nr*Nz),Nz,Nr);
vres = squeeze(vres(Nz/2,:)).'*exp(-1i*m*asangles);
ures = reshape(vv(Nr*Nz+1:end),Nz,Nr);
ures = squeeze(ures(Nz/2,:))'*exp(1i*m*asangles);
prres = squeeze(phi(N/2,N/2+1:end,Nz/2)).'*exp(1i*0*asangles) + 12*(vres + ures);
cmb = griddata(xas',yas',abs(prres).^2,grid.mesh.x2,grid.mesh.y2);
imagesc(grid.x*rscale*1e6,grid.y*rscale*1e6,cmb);
nv = sum(vdens*rr.')*(grid.x(2)-grid.x(1))*(grid.z(2)-grid.z(1))*2*pi;
% plot(rr,vdens(Nz/2,:)/nv);
drawnow;
% ddd=ddd(ddd>0);
end
% part1 = reshape(V+2*g*abs(phi).^2,N*N*Nz,1);
% part2 = reshape(g*phi.^2,N*N*Nz,1);
% K = 0.5*kron(kron(findiff1_arb(r),findiff1_arb(r)),findiff1_arb(rz));