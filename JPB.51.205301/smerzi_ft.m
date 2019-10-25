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
Ub=3500/omegascale*2*pi; %*1.5;
ub1=0;
NN = 1e5;
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
for k=1:1
    for i=1:15
        for j = 3:10
T = Tarr(j);
Tdim = kb*T/hbar/omegascale;

task1 = ZNGtask(grid,Vfun);
task1.g = g;
task1.Ntotal = NN;
task1.T = Tdim;
tshift = 0;
task1.Vtd = @TDPot_;

[phizng, mu, mu2] = task1.groundstate_itp(1e-2,1e-5,'tf');

mutot = gather(mu2(end)); %15.54*omegabar;

% [Nc,mutot1,ecut] = sgpe_params(omegabar,aRB/rscale,Tdim,8e4);
[Nc,~,ecut] = sgpe_params_mu(omegabar,aRB/rscale,Tdim,mutot);
ecut = real(ecut);
% mutot = real(mutot);
% RTF = sqrt(real(moutot)*2); %Thomas-Fermi radius
Veff = ecut-task1.getVtotal(0);

% ogrid = ogrid3d(1,1,omz,ecut,1);
ogrid2 = ogrid3d(1,1,omz,ecut,2);
Vfunp = @(X,Y,Z) X*0;
%%
task = PGPEtask(ogrid2,Vfunp);
task.g = 4*pi*aRB/rscale;
task.Ntotal = gather(grid.integrate(abs(phizng).^2));
task.n_crank = 5;
task.Vtd = @TDPot_;
ub1=k/10;
tshift = 1e10;
task.user_callback = @post_process_p;
task.gamma = real(4*mRB*aRB^2*kb*T/pi/hbar^2 *gamma_coef(mutot,ecut,Tdim));
[phitf, mutf] = task.groundstate_tf(1e-5);
% phitf(1:end/2,:,:) = -phitf(1:end/2,:,:);
[phi, mu, mu1] = task.groundstate_split(5e-3,1e-5,ogrid2.grid2sp(phitf));
phisp = task.init_state;
% task.ecut = mu2(end)*3;
%%
% j=0;
% Tarr = [278e-9, 250e-9, 300e-9, 350e-9, 400e-9, 450e-9];

% for gam = flip([0, 0.0002, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02, 0.03])
% for j = 1:1
%     T = Tarr(j);
    task.mu_init = mutot;
    task.n_recalc = 5;
    task.n_crank = 4;
    task.T = Tdim;
%     task.gamma = real(4*mRB*aRB^2*kb*T/pi/hbar^2 *gamma_coef(mutot,ecut,kb*T/hbar/omegascale));
%     ravg = 0;
%     j=j+1;

            task.current_iter = 0;
            task.snapshots={};
            task.init_state = phisp;
            tshift = 1e10;
            tic
            task.solve_crank(0.01,70,1000);
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
            save(sprintf('tasks_coh/task_%03d_%03d_%03d',j,k,i),'task');
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
times = (1:4000)*0.1*tcoef;%(1:4000)*0.04*tcoef;
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
for jj = 3:3
    ind=0;
    tmp = [];
for k=1:11

    G=0;
    phi2avg = 0;
    cnt = 0;
    for i=1:1%length(tss)
% %     if(exist(sprintf('tasks/task_%05d.mat',jj),'file') ~= 2)
% %         continue
% %     end
% %     load(sprintf('tasks/task_%05d.mat',jj));
    if(exist(sprintf('tasks/task_%03d_%03d_%03d.mat',jj,i,k),'file') ~= 2)
        continue
    end
    load(sprintf('tasks/task_%03d_%03d_%03d.mat',jj,i,k));
    ind=ind+1;
%     rr=sqrt(task.history.core1x.^2+task.history.core1y.^2);
%     ravg = ravg + rr;
%     cnt = cnt+1;
%     rr = ravg'/i;
%     tshift = tss(i)*1e-3;
%     plot(times,1-bar_pos(times,1,tshift)); drawnow;
% figure; hold all
    zz = (task.history.NN1-task.history.NN2)./(task.history.NN1+task.history.NN2);
    x = [ones(nfit(jj),1) times(sfit(jj)+1:nfit(jj)+sfit(jj))'];
    y = log(abs(zz(sfit(jj)+1:nfit(jj)+sfit(jj))'));
    b = x\y;
%     plot(times,-task.gamma*times*7/tcoef-0.5);
    
%     tmp(ind)=b(2);
    rr(jj)=(rr(jj)+b(2));
%     Nc(jj) = Nc(jj) + task.history.init.nc1(end);
%     Nc2(jj) = Nc2(jj) + task.history.init.nc2(end);
    gams1(jj)=task.gamma;
    ecuts(jj)=task.grid.ecut;
%     plot([Tarr(jj)/Tc],[gams1(jj)],'+');
%     plot([Tarr(jj)/Tc],[-b(2)*tcoef/lambdader(jj+10)],'x','Color',[0 0 0]);drawnow
    tmp1(end+1) = Tarr(jj)/Tc;
    tmp2(end+1) = -b(2)*tcoef/(Ua(jj+10)/1);
%     tmp2(end+1) = task.history.init.nc1(end);
    tmp3(end+1) = task.history.init.nc2(end);
    
%     plot(times, exp(b(1)+b(2)*times),'--','Color',[0 0 0]);
%     [zztm,phtm] = solve_two_mode(Ja(jj+10)/tcoef,0,Ua(jj+10)/1/tcoef,1/tmp2(end),exp(b(1)),0,times);
%     plot([times(sfit(jj)+1), times(sfit(jj)+1)],[-0.2,0.6],'--','Color',[0 0 0]);
%     plot([times(sfit(jj)+nfit(jj)), times(sfit(jj)+nfit(jj))],[-0.2,0.6],'--','Color',[0 0 0]);
%     plot(times,((zztm)),'Color',[0.85 0.33 0.1],'LineWidth',0.5); 
%      plot(times,((zz)),'Color',[0 0.45 0.74],'LineWidth',1); drawnow;
     
     plot(times,((zz))); drawnow;
%      [Nc,Nnc,~] = sgpe_params_mu(omegabar,aRB/rscale,task.T,task.mu_init);
%     plot([task.T],[gamma_coef(mutot,task.grid.ecut,task.T)],'o')
%      plot([task.T],[Nc],'o');
%      plot([task.T/(kb*Tc/hbar/omegascale)],[task.history.init.nc1(end)],'+'); drawnow;
%      plot([task.T/(kb*Tc/hbar/omegascale)],[task.history.init.nc2(end)],'o'); drawnow;
%      Nc = NN-polylog_inc(3,exp(task.mu_init/task.T),task.grid.ecut/task.T)/(omegabar/task.T)^3;
%      plot([task.T],[Nc],'x');
%      plot([task.T],[NN-Nnc],'s');
%     plot(task.history.init.nc1,'+')
%     plot(task.history.init.nc2,'.')
%     plot(task.history.nc1+task.history.nc2,'s')

%     phi4avg = 0;
%     for ii=1:length(task.snapshots)
%         G = G+task.snapshots{ii}*task.snapshots{ii}';
%         ssc = task.grid.etot*0;
%         ssc(task.grid.mask) = task.snapshots{ii};
%         sscr = task.grid.sp2grid(ssc);
%         phi2avg = phi2avg + abs(sscr).^2; 
%         cnt=cnt+1;
%     end
    end
%     phi2avg = phi2avg./cnt;
%     [current_c, current_nc] = eigs(G./cnt,10);
%     current_nc = diag(current_nc);
%     plot(current_nc); drawnow;
%     VVV = task.getVtotal(0);
%     MU2 = [];
%     phi = task.grid.etot*0;
%     phi(task.grid.mask) = current_c(:,1);
%     phir0 = task.grid.sp2grid(phi);   
%     phi = task.grid.etot*0;
%     phi(task.grid.mask) = current_c(:,2);
%     phir1 = task.grid.sp2grid(phi);     
%     for ii=1:10
%         phi = task.grid.etot*0;
%         phi(task.grid.mask) = current_c(:,ii);
%         phir = task.grid.sp2grid(phi);
%         MU2(ii) = real(sum(task.grid.to1d(task.grid.etot.*abs(phi).^2 + conj(phi).*(task.grid.grid2sp((VVV+task.g*(abs(phir0).^2*current_nc(1))).*phir)))));
%     end
%     plot(MU2);drawnow;
%     Nc(jj) = current_nc(1);
%     Nc2(jj) = current_nc(2);
    
%     ravg = ravg/cnt;
% 
% %     rr=sqrt(task.history.core1x.^2+task.history.core1y.^2)';
%     
%     theta=unwrap(atan2(task.history.core1x,task.history.core1y));
%     times=(1:numel(rr))'*2*tcoef;
% 
%     rr = rr(rr<13);
%     theta = theta(rr<13);
%     times = times(rr<13);   
%     
%     rr = rr(times>10);
%     theta = theta(times>10);
%     times = times(times>10);       
%     
%     rr = rr(1:step:end);
%     theta = theta(1:step:end);
%     times = times(1:step:end);
%     freq = findiff1_arb(times',1)*gather(theta');
%     
%     x = [ones(numel(rr),1) times];
%     y = log(rr);
%     bb = x\y;
%     plot(times,rr*rscale);
%     plot(times,exp(bb(1)+bb(2)*times)*rscale);
%     res(jj,1) = task.gamma;
%     res(jj,2) = gather(bb(2)/om0);
% %     plot(times(3:end-2),freq(3:end-2));
% %     plot(times(3:end-2),om0./(1-rr(3:end-2).^2/RTF^2)-1,'--');
% %     plot(task.history.core1x);
% %     plot(rr(3:end-2),freq(3:end-2))
%     
end
ind
rr(jj)=rr(jj)/ind;
sigma(jj) = sqrt(sum((tmp-rr(jj)).^2)/(ind-1));
end
% plot(Tarr/Tc,-rr*tcoef./lambdader(11:end-1),'s-')
% plot(Tarr/Tc,gams1,'s-')
% rr=rr/3;
% % plot(times,5*cos(om0./(1-rr.^2/RTF^2).*times/tcoef));
%%
Nctf = [];
Ncz = [];
Nncz = [];
Ncbarr = [];
muarr = [];
muarrz = [];
ecarr = [];
gams = [];
gams2 = [];
Tcdim = kb*Tc/omegascale/hbar;
for i=1:20
    task1 = ZNGtask(grid,Vfun);
    task1.g = g;
    task1.Ntotal = 1e5;
%     task1.mu_init = 10.407610298720511;
    tshift = 0;%1e10;
    task1.Vtd = @TDPot_;
    task1.T = i/20*Tcdim;
    [phizng, mu, mu2] = task1.groundstate_itp(1e-2,1e-5,'tf');
    Ncz(i) = gather(grid.integrate(abs(phizng).^2));
    Nncz(i) = gather(grid.integrate(abs(task1.init_state_nt)));
    muarrz(i) = gather(mu2(end));
    [Nctf(i),muarr(i),ecarr(i)] = sgpe_params(omegabar,aRB/rscale,task1.T,NN);
    [Ncbar(i),~,ecarr(i)] = sgpe_params_mu(omegabar,aRB/rscale,task1.T,muarrz(i));
    gams(i) = real(4*aRB^2*task1.T/pi/rscale^2 *gamma_coef(muarrz(i),ecarr(i),task1.T));
    gams2(i) = 1/task1.T*Nncz(i)/NN/2/pi*exp((muarrz(i)-Ub)/task1.T);
end

%%
gams2 = [];
gams3 = [];
ecarr = [];
Tcdim = kb*Tc/omegascale/hbar;
for i=1:20
    T=(i-1)/20*Tcdim;
    [~,~,ecarr(i)] = sgpe_params_mu(omegabar,aRB/rscale,T,MUtot(i));
    gams2(i) = 1/T*tmp4(i)/NN/2/pi*exp((MUtot(i)-Ub)/T);
    gams3(i) = real(4*aRB^2*T/pi/rscale^2 *gamma_coef(MUtot(i),ecarr(i),T));
end

%%



figure;

times = (1:4000)*0.1*tcoef;%(1:4000)*0.04*tcoef;
ravg = 0;
cnt = 0;
rr=zeros(1,9);
sigma=zeros(1,9);
tmp = [];
tmp1 = [];
tmp2 = [];
tmp3 = [];
gams1=[];
ecuts =[];
nfit = [0, 0, 900, 1100, 1400, 1400, 2100, 3500, 3400];
sfit = [0, 0, 600, 700, 1100, 1100, 1300, 300, 400];
ax = [];
ax(1) = axes('Parent',gcf,'FontSize',9,...
    'Position',[0.08 0.1 0.4 0.38]);
ax(2) = axes('Parent',gcf,'FontSize',9,...
    'Position',[0.08 0.6 0.4 0.38]);%,'XTickLabel',[]);
ax(3) = axes('Parent',gcf,'FontSize',9,...
    'Position',[0.58 0.1 0.4 0.38]);%,'YTickLabel',[]);
ax(4) = axes('Parent',gcf,'FontSize',9,...
    'Position',[0.58 0.6 0.4 0.38]);%,'YTickLabel',[],'XTickLabel',[]);
xlim(ax(1),[0 0.28]);
xlim(ax(2),[0 0.28]);
xlim(ax(3),[0 0.28]);
xlim(ax(4),[0 0.28]);
ylim(ax(1),[-0.1 0.6]);
ylim(ax(2),[-0.1 0.6]);
ylim(ax(3),[0 0.6]);
ylim(ax(4),[-0.1 0.6]);

ylabel(ax(1),'$Z$','FontSize',12,'Interpreter','latex');
ylabel(ax(2),'$Z$','FontSize',12,'Interpreter','latex');
ylabel(ax(3),'$Z$','FontSize',12,'Interpreter','latex');
ylabel(ax(4),'$Z$','FontSize',12,'Interpreter','latex');

xlabel(ax(1),'time (s)','FontSize',12,'Interpreter','latex');
xlabel(ax(3),'time (s)','FontSize',12,'Interpreter','latex');
xlabel(ax(2),'time (s)','FontSize',12,'Interpreter','latex');
xlabel(ax(4),'time (s)','FontSize',12,'Interpreter','latex');

box(ax(1),'on');
box(ax(2),'on');
box(ax(3),'on');
box(ax(4),'on');
ind=1;
shows = [2,8,9,5];
for jj = [7,3,9,5]
    
    tmp = [];
    set(gcf, 'currentaxes', ax(ind));
    hold on;
for k=[(1:11) shows(ind)]

    G=0;
    cnt = 0;
    for i=1:1%length(tss)

    if(exist(sprintf('tasks/task_%03d_%03d_%03d.mat',jj,i,k),'file') ~= 2)
        continue
    end
    load(sprintf('tasks/task_%03d_%03d_%03d.mat',jj,i,k));
    

% figure; hold all
    zz = (task.history.NN1-task.history.NN2)./(task.history.NN1+task.history.NN2);
    x = [ones(nfit(jj),1) times(sfit(jj)+1:nfit(jj)+sfit(jj))'];
    y = log(abs(zz(sfit(jj)+1:nfit(jj)+sfit(jj))'));
    b = x\y;

    rr(jj)=(rr(jj)+b(2));

    gams1(jj)=task.gamma;
    ecuts(jj)=task.grid.ecut;

    tmp1(end+1) = Tarr(jj)/Tc;
    tmp2(end+1) = -b(2)*tcoef/(Ua(jj+10)/1);

%     if(jj==3 || jj==7)
    plot(times,((zz)),'Color',[0.8 0.8 0.8],'LineWidth',0.5); drawnow;
%     end
    end

end
%     plot(times, exp(b(1)+b(2)*times),'--','Color',[0 0 0]);
    [zztm,phtm] = solve_two_mode(Ja(jj+10)/tcoef,0,Ua(jj+10)/1/tcoef,1/tmp2(end),exp(b(1)),0,times);
    plot([times(sfit(jj)+1), times(sfit(jj)+1)],[-0.2,0.6],'--','Color',[0 0 0]);
    plot([times(sfit(jj)+nfit(jj)), times(sfit(jj)+nfit(jj))],[-0.2,0.6],'--','Color',[0 0 0]);
 plot(times(1:end),((zztm(1:end))),'Color','blue','LineWidth',0.5);
     plot(times,((zz)),'Color','red','LineWidth',0.5); drawnow;
          
     ind=ind+1;
end

annotation(gcf,'textbox',...
    [0.3 0.9 0.16 0.06],...
    'String',{'$T=0.6T_\mathrm{c}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FontName','Helvetica');
annotation(gcf,'textbox',...
    [0.85 0.9 0.16 0.06],...
    'String',{'$T=0.7T_\mathrm{c}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FontName','Helvetica');
annotation(gcf,'textbox',...
    [0.27 0.4 0.16 0.06],...
    'String',{'$T=0.8T_\mathrm{c}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FontName','Helvetica');
annotation(gcf,'textbox',...
    [0.817 0.4 0.16 0.06],...
    'String',{'$T=0.9T_\mathrm{c}$'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',12,...
    'FontName','Helvetica');
