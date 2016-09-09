%% calculate coefficients of TM approx


global tcoef Ub omegaBar initphase ub1;
%   Box size
xmax = 15;
zmax = 7.5;
%   Grid size
N=128;
Nz=64;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);

tcoef = 1/(570*2*pi); %time scale
Ub=1.60931; %*1.5
ub1 = 0;
omegaBar = 0.0;
initphase = 0 ;%-0.125*pi;

r0 = 8.85787;
omz = (0.526316);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + (sqrt(X.^2+Y.^2)-r0).^2);

task = GPEtask(grid,Vfun);
task.g = 0.160567;%/sqrt(2*pi/omz);
task.Ntotal = 2000; 
task.Vtd = @TDPot;

Ja = zeros(1,30);
Jb = zeros(1,30);
K2 = zeros(1,30);
FF = zeros(1,30);
Ua = zeros(1,30);
Udif = zeros(1,30);
coeffs = zeros(1,30);
E1 = zeros(1,30);
E2 = zeros(1,30);
MUtot = zeros(1,30);
nd = numel(size(grid.mesh.x));
if(nd==3)
    weightxz = (grid.x(2)-grid.x(1))*(grid.z(2)-grid.z(1));
else
    weightxz = (grid.x(2)-grid.x(1));
end
angles3 = atan2(grid.mesh.x,grid.mesh.y);
%%
Ntotal = task.Ntotal;
for jj = 1:1
% load('phi_1000_1.mat');
% task.Ntotal = 250*jj;
% initphase = 0.02*(jj-1);
% deltaN = (jj-1)*Ntotal/30;
% NN0 = Ntotal + deltaN;
% Ub = jj*0.1;
[phiplus,mu1] = task.groundstate_itp(1e-2,1e-6);
if(nd==3)
	coeffs(jj)=gather(sum(squeeze(abs(phiplus(N/2,N/2:end,:)).^2)'*grid.x(N/2:end)')*weightxz/task.Ntotal*4);
else
    coeffs(jj)=gather(abs(phiplus(N/2,N/2:end)).^2*grid.x(N/2:end)'*weightxz/task.Ntotal*4);
end
MUtot(jj) = gather(mu1(end));
% load('phi_1000_2.mat');
phi0 = phiplus;
% phi0(:,1:N/2,:)=-phi0(:,1:N/2,:);
phi0(angles3>(-pi + initphase) & angles3<(- initphase)) = -phi0(angles3>(-pi + initphase) & angles3<(- initphase));
[phimin,mu2] = task.groundstate_itp(1e-2,1e-6,phi0);

% NN0 = Ntotal - deltaN;
% itp_run;
% phiplus2=phi;
% coeffs(jj)=abs(sum(squeeze(abs(phi(64,64:end,:)).^2)'*r(64:end)')*h*hz/NN0*4);
% % load('phi_1000_2.mat');
% phi0 = phi;
% % phi0(:,1:N/2,:)=-phi0(:,1:N/2,:);
% phi0(angles3>(-pi + initphase) & angles3<(- initphase)) = -phi0(angles3>(-pi + initphase) & angles3<(- initphase));
% itp_run;
% phimin2=phi;


% NN0=sum(sum(sum(abs(phi).^2)))*h*h*hz
phi1 = (phiplus + phimin)/sqrt(2);
phi2 = (phiplus - phimin)/sqrt(2);
N1 = grid.integrate(abs(phi1).^2);
N2 = grid.integrate(abs(phi2).^2);
phi1 = phi1/sqrt(N1);
phi2 = phi2/sqrt(N2);
% imagesc(real(phi2(:,:,64))); drawnow;

V = task.getVtotal(0);
K = -grid.inner(phi1,grid.ifft(grid.kk.*grid.fft(phi2)) + V.*phi2);
F12 = -task.g*grid.inner(phi1,abs(phi1).^2.*phi2);
F21 = -task.g*grid.inner(phi1,abs(phi2).^2.*phi2);
U1 = task.g*grid.inner(phi1,abs(phi1).^2.*phi1);
U2 = task.g*grid.inner(phi2,abs(phi2).^2.*phi2);

II = task.g*grid.integrate(abs(phi2).^2.*abs(phi1).^2);
E1 = gather(grid.inner(phiplus,grid.ifft(grid.kk.*grid.fft(phiplus)) + V.*phiplus + task.g*abs(phiplus).^2.*phiplus))
E2 = gather(grid.inner(phimin,grid.ifft(grid.kk.*grid.fft(phimin)) + V.*phimin + task.g*abs(phimin).^2.*phimin))

K2(jj) = gather(2*K);
FF(jj) = gather(F12+F21);
Ja(jj) = gather(real(2*K+task.Ntotal*(F12+F21)));
Jb(jj) = gather(real(task.Ntotal*(F21-F12)));
Ua(jj) = gather(real((U1+U2)/2));
Udif(jj) = gather(real((U1-U2)/2));
%
% times = (0:5000)*0.5;
% zz0=coeffs(jj)*(arrayfun(@rotangle,times,0.5+times*0));
% [zz,dphi] = solve_two_mode(Ja(jj),Jb(jj),Ua(jj)*task.Ntotal,400,0,0,times,zz0);

end