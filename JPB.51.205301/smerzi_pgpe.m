% global tcoef Ub omegaBar initphase ub1 r0 r1 r01;
units
global linit tcoef Ub rscale vv ub1 tshift;


% Smerzi PhysRevLett.84.4521
omegascale = 50*2*pi;
omz = 17.68/50;
omegabar = omz^(1/3);
tcoef = 1/(omegascale); %time scale
rscale = sqrt(hbar/mRB/omegascale);
linit=0;
vv = 0;%(4e-7)/(rscale/tcoef);
Ub=650/omegascale*2*pi; %*1.5;
ub1=0;
NN = 5.0e4;
%   Box size
xmax = 20;
zmax = 20;
%   Grid size
N=64;
Nz=64;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);

g = 4*pi*aRB/rscale;% /sqrt(2*pi/omz);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + X.^2 + Y.^2);

%%

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
NN = 5e4;
%   Box size
xmax = 17;
zmax = 90;
%   Grid size
N=64;
Nz=128;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);

g = 4*pi*aRB/rscale;% /sqrt(2*pi/omz);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + X.^2 + Y.^2);

%%

%     task1 = GPEtask(grid,Vfun);
%     task1.g = g;
%     task1.Ntotal = NN;
%     task1.Vtd = @TDPot;
%     [phi, mu, mu2] = task1.groundstate_itp(1e-2,1e-5,'tf');
    
%%

% T = 230e-9; %in K
% 
% task1 = SGPEtask(grid,Vfun);
% task1.g = g;
% task1.Ntotal = 8.0e4; 
% % task.mu_init = 10;
% % task.Vtd = @TDPot_c;
% % task.gamma = 0.03;
% task1.user_callback = @post_process;
% % task.T = 85e-9 *kb/hbar/omegascale;
% % task.ecut = task.T*log(4/3) + task.mu_init;
% 
% [phi1, mu3] = task1.groundstate_tf(1e-6);
% [phi, mu, mu2] = task1.groundstate_itp(1e-3,1e-6,phi1);
% xi = gather(1/sqrt(2*g*max(abs(phi(:))).^2)); % healing length
% RTF = gather(sqrt(real(mu2(end))*2)); %Thomas-Fermi radius
% 
% task1.mu_init = mu2(end);
% task1.n_recalc=1000;
% mutot = real(gather(mu2(end)));

%%
Tc=225e-9; % Levi
% figure;
Tarr = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]*Tc;
tss = [5,10,15,20,25,30,35,40,45,50,60,80,100];
%% initial parameter estimate for the fixed-mu case
% T = Tarr(1);
% Tdim = kb*T/hbar/omegascale;
% task1 = ZNGtask(grid,Vfun);
% task1.g = g;
% task1.Ntotal = NN;
% task1.T = Tdim;
% tshift = 0;
% task1.Vtd = @TDPot_;
% 
% [phizng, mu, mu2] = task1.groundstate_itp(1e-2,1e-5,'tf');
% 
% mutot = gather(mu2(end)); %15.54*omegabar;
% ecut = Tdim*log(1+1/2) + mutot; %thermal population of 2 atoms per mode at Ecut %26.35*omegabar;
% Nnc = polylog_inc(3,exp(mutot/Tdim),ecut/Tdim)/(omegabar/Tdim)^3;
% xi = gather(1/sqrt(2*g*max(abs(phizng(:))).^2)); % healing length
%% calculating the density of states in full HF approximation
% U = task1.getVtotal(0)+2*task1.g*(abs(phizng).^2 + abs(task1.init_state_nt))-mutot;
% % U = abs(task1.getVtotal(0)-mutot)+mutot;
% rho1 = zeros(1,1000);
% rho1int = zeros(1,1000);
% ecut1 = zeros(1,1000);
% ecut2 = zeros(1,1000);
% Tcdim = kb*Tc/hbar/omegascale;
% emax = Tcdim*log(2);
% de = emax/1000;
% e0=1+omz/2;
% for i=1:1000
%     e=i*de;
%     TT = Tcdim*i/1000;
% %     mask = U>e;
% %     U1 = U;
% %     U1(mask) = 0;
%     rho1(i)= gather(grid.integrate(real(sqrt(e-gather(U))))/sqrt(2)/pi^2);
%     rho1int(i) = sum(rho1)*de;
%     ecut1(i) = (rho1int(i)*6*omegabar^3 + e0^3)^(1/3);
%     [~,~,ecut2(i)] = sgpe_params_mu(omegabar,aRB/rscale,TT,mutot);
% end
%%
for k=1:1%length(tss)
    for i=2:2
        for j = 1:10
% T = Tarr(j);
% Tdim = kb*T/hbar/omegascale;
% 
% task1 = ZNGtask(grid,Vfun);
% task1.g = g;
% task1.Ntotal = NN;
% task1.T = Tdim;
% tshift = 0;
% task1.Vtd = @TDPot_;
% 
% [phizng, mu, mu2] = task1.groundstate_itp(1e-2,1e-5,'tf');
% 
% mutot = gather(mu2(end)); %15.54*omegabar;

% [Nc,mutot1,ecut] = sgpe_params(omegabar,aRB/rscale,Tdim,8e4);
% [Nc,~,ecut] = sgpe_params_mu(omegabar,aRB/rscale,Tdim,mutot);
ecut = j+9;
% mutot = real(mutot);
% RTF = sqrt(real(moutot)*2); %Thomas-Fermi radius
% Veff = ecut-task1.getVtotal(0);

% ogrid = ogrid3d(1,1,omz,ecut,1);
ogrid2 = ogrid3d(1,1,omz,ecut,1);
Vfunp = @(X,Y,Z) X*0;
%%
task = PGPEtask(ogrid2,Vfunp);
task.g = 4*pi*aRB/rscale;
task.Ntotal = NN;
task.n_crank = 5;
task.Vtd = @TDPot_;
% ub1=k/10;
tshift = 0;%1e10;
task.user_callback = @post_process_p;
task.gamma = 0;%real(4*mRB*aRB^2*kb*T/pi/hbar^2 *gamma_coef(mutot,ecut,Tdim));
[phitf, mutf] = task.groundstate_tf(1e-5);
% phitf(1:end/2,:,:) = -phitf(1:end/2,:,:);
[phi, mu, mu1] = task.groundstate_split(5e-3,1e-6,ogrid2.grid2sp(phitf));
phisp = task.init_state;
% task.ecut = mu2(end)*3;
%%
% j=0;
% Tarr = [278e-9, 250e-9, 300e-9, 350e-9, 400e-9, 450e-9];

% for gam = flip([0, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03])
% for j = 1:1
%     T = Tarr(j);
    task.mu_init = 0;
    task.n_recalc = 5;
    task.n_crank = 4;
    task.T = 0;
%     task.gamma = real(4*mRB*aRB^2*kb*T/pi/hbar^2 *gamma_coef(mutot,ecut,kb*T/hbar/omegascale));
%     ravg = 0;
%     j=j+1;

            task.current_iter = 0;
            task.snapshots={};
            task.init_state = phisp;
%             tshift = tss(k)*1e-3/tcoef;
            tic
            task.solve_split(0.01,50,1000);
%             task.history.init = task.history;
            toc  
    %         phi2 = phisp.*0;
    %         phi2(ogrid2.mask) = task.current_c;
    %         phi2r = ogrid2.sp2grid(phi2)*sqrt(task.current_nc);
            %%
%             task.current_iter = 0;
%             task.init_state = task.current_state;
    %         phiv = ogrid2.sp2grid(task.current_state).*tanh(sqrt((ogrid2.mesh.x-0.5).^2 + ogrid2.mesh.y.^2)./xi).^1.*exp(1i*atan2(ogrid2.mesh.x-0.5,ogrid2.mesh.y));
    %         task.init_state = ogrid2.grid2sp(phiv);
    % % %         task.init_state = task.current_state.*tanh(sqrt((grid.mesh.x-1).^2 + grid.mesh.y.^2)./xi).^1.*exp(1i*atan2(grid.mesh.x-1,grid.mesh.y));
    %         tic
%             tshift = tss(1)*1e-3/tcoef;
%             task.solve_crank(0.005,20,4000);
    %         toc
    %         rr=sqrt(task.history.core1x.^2+task.history.core1y.^2);
    %         ravg = ravg + rr;
            save(sprintf('tasks_silent/task_%03d_%03d_%03d',j,k,i),'task');
        end
    end
end

%%

%     G=0;
%     for i=1:length(task.snapshots)
%         G = G+task.snapshots{i}*task.snapshots{i}';
%     end
%     %%
%     [V, NNC] = eigs(G./length(task.snapshots),1);

% %%
figure; hold all;
% om0=3/2/RTF^2*log(RTF/xi)/tcoef;
% res = zeros(6,2);
% step = 5;
times = (1:1000)*0.5*tcoef;%(1:4000)*0.04*tcoef;
ravg = 0;
cnt = 0;
% Nc = zeros(1,9);
% Nc2 = zeros(1,9);
rr=zeros(1,9);
sigma=zeros(1,9);
tmp = [];
tmp1 = [];
tmp2 = [];
tmp3 = [];
gams1=[];
ecuts =[];
nfit = [0, 0, 900, 1100, 1400, 1400, 2000, 3500, 3500];
sfit = [0, 0, 600, 700, 1100, 1100, 1400, 300, 300];
for jj = 1:10
    ind=0;
    tmp = [];
for k=2:2

    G=0;
    phi2avg = 0;
    cnt = 0;
    for i=1:10
% %     if(exist(sprintf('tasks/task_%05d.mat',jj),'file') ~= 2)
% %         continue
% %     end
% %     load(sprintf('tasks/task_%05d.mat',jj));
    if(exist(sprintf('tasks_silent/task_%03d_%03d_%03d.mat',jj,i,k),'file') ~= 2)
        continue
    end
    load(sprintf('tasks_silent/task_%03d_%03d_%03d.mat',jj,i,k));
    ind=ind+1;
%     rr=sqrt(task.history.core1x.^2+task.history.core1y.^2);
%     ravg = ravg + rr;
%     cnt = cnt+1;
%     rr = ravg'/i;
%     tshift = tss(i)*1e-3;
%     plot(times,1-bar_pos(times,1,tshift)); drawnow;
% figure; hold all
    zz = (task.history.NN1-task.history.NN2)./(task.history.NN1+task.history.NN2);

     plot(times,((task.history.n))); drawnow;

    end

end

end

