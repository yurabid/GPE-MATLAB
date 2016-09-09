global tcoef Ub omegaBar initphase ub1;
%   Box size
xmax = 15;
zmax = 7.5;
%   Grid size
N=128;
Nz=32;
grid = grid2d(xmax,N,xmax,N,zmax,Nz);

tcoef = 1/(570*2*pi); %time scale
Ub=1.60931; %*1.5
ub1 = 0.0;
omegaBar = 0.0;
initphase = 0.0;%-0.06*pi;



r0 = 8.85787;
omz = (0.526316);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + (sqrt(X.^2+Y.^2)-r0).^2);

task = GPEtask(grid,Vfun);
task.g = 0.160567/sqrt(2*pi/omz);
task.Ntotal = 5000; 
% task.Vtd = @TDPot;
task.omega = 0.0;
s = 0;
[phi, mu1, mu2] = task.groundstate_itp(1e-3,1e-7);
% phi = phi.*exp(-1i*s*atan2(grid.mesh.x,grid.mesh.y));
% [phi, mu1, mu2] = task.groundstate_itp(2e-3,1e-7,phi);
% phi = phi.*exp(-1i*15*atan2(grid.mesh.x,grid.mesh.y));
mu=grid.inner(phi,task.applyham(phi))/task.Ntotal; %mu1(end);
eps = task.applyham(phi) - mu*phi;
%%

spectr = zeros(60,20);
% figure;
hold all;
for i=1:60
m=i-1;
rr=grid.x(N/2+1:end);
Nr = length(rr);
% [RR,ZZ] = meshgrid(rr,grid.z);
Vr = reshape(Vfun(rr,0,0),Nr,1);
phir = abs(squeeze(phi(N/2,N/2+1:end)));
lapr = findiff1_arb(rr,2,0) + spdiags(1./rr(:),0,Nr,Nr)*findiff1_arb(rr,1,0);

nlin = task.g*spdiags(abs(phir(:)).^2,0,Nr,Nr);
h1 = -0.5*(lapr - (m+s)^2*spdiags(1./rr(:).^2,0,Nr,Nr)) + spdiags(Vr(:),0,Nr,Nr) + 2*nlin - mu*speye(Nr);
h2 = -0.5*(lapr - (m-s)^2*spdiags(1./rr(:).^2,0,Nr,Nr)) + spdiags(Vr(:),0,Nr,Nr) + 2*nlin - mu*speye(Nr);
mtot = [h1 nlin; -nlin -h2];
[vv,dd]=eigs(mtot,1,'sm');
ddd=diag(dd);
spectr(i,:) = ddd(:);
plot(m,(abs(ddd)),'LineStyle','none','Marker','.','Color',[1 0 0]);
% plot((1:20)*0+m,(imag(ddd)),'LineStyle','none','Marker','.','Color',[1 0 0]);
drawnow;
% ddd=ddd(ddd>0);
end
% part1 = reshape(V+2*g*abs(phi).^2,N*N*Nz,1);
% part2 = reshape(g*phi.^2,N*N*Nz,1);
% K = 0.5*kron(kron(findiff1_arb(r),findiff1_arb(r)),findiff1_arb(rz));