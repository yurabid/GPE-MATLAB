%% calculate coefficients of TM approx

units
global linit tcoef Ub rscale vv ub1 tshift;

% Levi
omegascale = 224*2*pi;
omz = 26/224;
omegabar = omz^(1/3);
tcoef = 1/(omegascale); %time scale
rscale = sqrt(hbar/mRB/omegascale);
linit=0;
vv = 0;%(4e-7)/(rscale/tcoef);
Ub=3000/omegascale*2*pi; %*1.5;
ub1=0;
NN = 1e5;
%   Box size
xmax = 17;
zmax = 90;
%   Grid size
N=64;
Nz=128;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);
grid2 = grid3dgpu(10,N,10,N,50,Nz);
%%
g = 4*pi*aRB/rscale;% /sqrt(2*pi/omz);
% Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2  + Y.^2) + Ub*exp(-2*((X )./(0.7e-6/rscale)).^2).*exp(-2*((Y )./(200e-6/rscale)).^2);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + X.^2 + Y.^2);
tshift = 0;
Tc=225e-9 *kb/omegascale/hbar; % Levi
Tarr = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]*Tc;

Ja = zeros(1,9);
Jb = zeros(1,9);
K2 = zeros(1,9);
FF = zeros(1,9);
Ua = zeros(1,9);
Ua2 = zeros(1,9);
Udif = zeros(1,9);
coeffs = zeros(1,9);
E1 = zeros(1,9);
E2 = zeros(1,9);
Nc = zeros(1,9);
Nc1 = zeros(1,9);
MUtot = zeros(1,9);
lambdader = zeros(1,9);
tmp = zeros(1,9);
tmp2 = zeros(1,9);
tmp3 = zeros(1,9);
tmp4 = zeros(1,9);
nd = numel(size(grid.mesh.x));
if(nd==3)
    weightxz = (grid.x(2)-grid.x(1))*(grid.z(2)-grid.z(1));
else
    weightxz = (grid.x(2)-grid.x(1));
end
angles3 = atan2(grid.mesh.x,grid.mesh.y);

%%
% mutf = [];
% ecuttf=[];
% Nctf=[];
% for i=1:200
% [Nctf(i),mutf(i),ecuttf(i)] = sgpe_params(omegabar,aRB/rscale,i/200*Tc,NN);
% end

%%
%Ntotal = task.Ntotal;
figure;
hold all


% dmu = 0.02;
for jj = 1:20
    
    task = ZNGtask(grid,Vfun);
    task.g = g;
    task.Vtd = @TDPot_;
    task.T=(jj-1)/20*Tc;
%     task.mu_init = 0;%10.4992-dmu;
    task.Ntotal = NN-1000;
    [phiplus1,mu11,mu22] = task.groundstate_itp(1e-2,1e-5);
    Nc1(jj) = gather(abs(grid.inner(phiplus1,phiplus1)));
    tmp3(jj) = gather(grid.integrate(task.init_state_nt));
    task.mu_init = 0;%10.4992;
    task.Ntotal = NN;
%     task.Ntotal = NN-jj*1000;
% load('phi_1000_1.mat');
% task.Ntotal = 250*jj;
% initphase = 0.02*(jj-1);
% deltaN = (jj-1)*Ntotal/30;
% NN0 = Ntotal + deltaN;
% Ub = jj*0.1;
[phiplus,mu1,mu2] = task.groundstate_itp(1e-2,1e-5,phiplus1);
ntplus = task.init_state_nt;
Nc(jj) = gather(abs(grid.inner(phiplus,phiplus)));
tmp4(jj) = gather(grid.integrate(task.init_state_nt));
lambdader(jj) = gather((Nc(jj)+Nc1(jj))*(mu22(end)-mu2(end))/(Nc1(jj)-Nc(jj)));
% [Nctf1,mutottf1,~] = sgpe_params(omegabar,aRB/rscale,task.T,NN-1000);
% [Nctf2,mutottf2,~] = sgpe_params(omegabar,aRB/rscale,task.T,NN);
% lambdader(jj) = gather(mu2(end));
% tmp(jj) = gather((Nctf1+Nctf2)*(mutottf1-mutottf2)/(Nctf1-Nctf2));
% tmp2(jj) = gather((Nc(jj)+Nc1(jj))*(mu22(end-1)-mu2(end-1))/(Nc1(jj)-Nc(jj)));
% plot(mu1(end-500:end));
% plot(real(mu2(end-500:end)));drawnow;
% task.Ntotal = Nc(jj);
% U = task.getVtotal(0)+2*task.g*(abs(phiplus).^2 + abs(task.init_state_nt))-mu2(end);
% rho1 = zeros(1,1000);
% rho1int = zeros(1,1000);
% ecut1 = zeros(1,1000);
% ecut2 = zeros(1,1000);
% Tcdim = kb*Tc/hbar/omegascale;

% emax = task.T*log(2);
% de = emax/1000;
% e0=1+omz/2;
% for i=1:1000
%     e=i*de;
%     rho1(i)= gather(grid.integrate(real(sqrt(e-gather(U))))/sqrt(2)/pi^2);
% end
%     rho1int = sum(rho1)*de;
%     tmp2(jj) = (rho1int*6*omegabar^3 + e0^3)^(1/3);
    
% if(nd==3)
% 	coeffs(jj)=gather(sum(squeeze(abs(phiplus(N/2,N/2:end,:)).^2)'*grid.x(N/2:end)')*weightxz/task.Ntotal*4);
% else
%     coeffs(jj)=gather(abs(phiplus(N/2,N/2:end)).^2*grid.x(N/2:end)'*weightxz/task.Ntotal*4);
% end
phiplus0 = phiplus/sqrt(Nc(jj));
phiplus1 = phiplus1/sqrt(Nc1(jj));
tmp(jj) = gather(grid.integrate(abs(phiplus0).^2.*abs(phiplus1).^2)/grid.integrate(abs(phiplus0).^4))-1;
tmp3(jj) = gather(grid.integrate(abs(phiplus1).^2.*abs(phiplus1).^2)/grid.integrate(abs(phiplus0).^4))-1;
MUtot(jj) = gather(mu2(end));
% tmp3(jj) = gather(mu1(end));
Ua(jj) = gather(grid.integrate(abs(phiplus0).^4)*task.g*Nc(jj)*2);
tmp2(jj) = Ua(jj).*(1-tmp(jj).*Nc(jj)./(Nc(jj)-Nc1(jj)));
jj
lambdader(jj)
tmp2(jj)
% tmp(jj)
% tmp(jj)/(Nc(jj)-Nc1(jj))*Nc(jj)
% tmp2(jj)
% MUtot(jj)
% tmp(jj)
% continue
% load('phi_1000_2.mat');

% task = GPEtask(grid2,Vfun);
% task.g = g;
% task.Ntotal = Nc(jj); 
% task.Vtd = @TDPot_;
% [phiplus,mu1] = task.groundstate_itp(1e-2,1e-5);
phiplus = abs(phiplus);
phi0 = phiplus;
phi0(:,1:N/2,:)=-phi0(:,1:N/2,:);
% phi0(angles3>(-pi + initphase) & angles3<(- initphase)) = -phi0(angles3>(-pi + initphase) & angles3<(- initphase));
[phimin,mu2] = task.groundstate_itp(1e-2,1e-5,phi0);
ntmin = task.init_state_nt;
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
N1 = grid2.integrate(abs(phi1).^2);
N2 = grid2.integrate(abs(phi2).^2);
phi1 = real(phi1)/sqrt(N1);
phi2 = real(phi2)/sqrt(N2);
% imagesc(real(phi2(:,:,64))); drawnow;

V = task.getVtotal(0);
K = -grid2.inner(phi1,grid2.ifft(grid2.kk.*grid2.fft(phi2)) + (V+2*g*ntplus).*phi2);
F12 = -task.g*grid2.inner(phi1,abs(phi1).^2.*phi2);
F21 = -task.g*grid2.inner(phi1,abs(phi2).^2.*phi2);
U1 = task.g*grid2.inner(phi1,abs(phi1).^2.*phi1);
U2 = task.g*grid2.inner(phi2,abs(phi2).^2.*phi2);

II = task.g*grid2.integrate(abs(phi2).^2.*abs(phi1).^2);
E1(jj) = gather(grid2.inner(phiplus,grid2.ifft(grid2.kk.*grid2.fft(phiplus)) + (V+2*g*ntplus).*phiplus + task.g/2*abs(phiplus).^2.*phiplus));
E2(jj) = gather(grid2.inner(phimin,grid2.ifft(grid2.kk.*grid2.fft(phimin)) + (V+2*g*ntmin).*phimin + task.g/2*abs(phimin).^2.*phimin));

Nc(jj) = gather(abs(grid2.inner(phiplus,phiplus)));
K2(jj) = gather(2*K);
FF(jj) = gather(F12+F21);
Ja(jj) = gather(real(2*K+Nc(jj)*(F12+F21)));
Jb(jj) = gather(real(task.Ntotal*(F21-F12)));
Ua2(jj) = gather(real((U1+U2)/2));
Udif(jj) = gather(real((U1-U2)/2));
%
% times = (0:5000)*0.5;
% zz0=coeffs(jj)*(arrayfun(@rotangle,times,0.5+times*0));
% [zz,dphi] = solve_two_mode(Ja(jj),Jb(jj),Ua(jj)*task.Ntotal,400,0,0,times,zz0);
Ua2(jj)*Nc(jj)*7/10
Ja(jj)*Nc(jj)/2

end
%%
% ecut2=[];
% for jj = 1:20
% [~,~,ecut2(jj)] = sgpe_params_mu(omegabar,aRB/rscale,(jj-1)/20*Tc,MUtot(jj));
% end

%%
% task = GPEtask(grid2,Vfun);
% task.g = g;
% task.gamma = 1.3e-3;
% tshift = 5*1e-3/tcoef;
% task.Ntotal = NN/2; 
% task.Vtd = @TDPot_;
% task.user_callback = @post_process;
% [phiplus,mu1] = task.groundstate_itp(1e-2,1e-5);
% task.solve_split(0.005,20,4000)