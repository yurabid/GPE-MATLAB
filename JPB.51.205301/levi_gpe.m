% global tcoef Ub omegaBar initphase ub1 r0 r1 r01;
units
global linit tcoef Ub rscale vv ub1 tshift;


%%
% amps = [4000];
% for j = 1:length(amps)
% Levi
omegascale = 224*2*pi;
omz = 26/224;
omegabar = omz^(1/3);
tcoef = 1/(omegascale); %time scale
rscale = sqrt(hbar/mRB/omegascale);
linit=0;
vv = 0;%(4e-7)/(rscale/tcoef);
Ub=3000/omegascale*2*pi; %*1.5;
ub1=0.3;
NN = 1e5;
%   Box size
xmax = 8;
zmax = 50;
%   Grid size
N=64;
Nz=128;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);

g = 4*pi*aRB/rscale;% /sqrt(2*pi/omz);
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + X.^2 + Y.^2);
% Vfun = @(X,Y,Z) 0.5*(0*omz^2*Z.^2  + Y.^2) + 20 - 20./(1+exp(-Z-40)) + 20./(1+exp(-Z+40)) + Ub*exp(-2*((X )./(0.7e-6/rscale)).^2);%.*exp(-2*((Z )./(200e-6/rscale)).^2);
tss = [5,10,15,20,25,30,35,40,45,50,60,80,100];
%%
    task1 = GPEtask(grid,Vfun);
    tshift = 10;%tss(2)*1e-3/tcoef;
    task1.g = g;
    task1.Ntotal = NN;
    task1.gamma=0.00;
%     task1.n_recalc=5;
    task1.Vtd = @TDPot2;
    task1.user_callback = @post_process;
    mus = [];
% for i=1:100
%     task1.Ntotal = NN*i/100;
    [phi, mu, mu2] = task1.groundstate_itp(1e-2,1e-5,'tf');
%     mus(i) = gather(mu2(end));
% end
%     ang = atan2(grid.mesh.x,grid.mesh.z-3);
%     phi = phi.*exp(1i*ang);
%     [phi, mu3, mu4] = task1.groundstate_itp(1e-2,1e-5,phi);
%     phi(:,1:end/2,:) = phi(:,1:end/2,:).*exp(2i*atan(exp(0.3*(grid.mesh.z(:,1:end/2,:)-3))));
%     phi(:,1:end/2,:) = phi(:,1:end/2,:).*atan(((grid.mesh.z(:,1:end/2,:)*10)))*2/pi;
%      phi(:,end/2+1:end,:) = 0.90*phi(:,end/2+1:end,:).*exp(-0i*atan(exp(0.3*(grid.mesh.z(:,end/2+1:end,:)-3))));
%     [phi, mu3, mu4] = task1.groundstate_itp(1e-2,1e-5,phi);
%     phi(:,1:end/2,:) = 0.95*phi(:,1:end/2,:);
%     task1.init_state = phi;
    %%
    for i=13:13%length(tss)
%         tshift = 0;%tss(i)*1e-3/tcoef;
        task1.current_iter=0;
        tic
        task1.solve_split(0.005,20,4000);
%         save(sprintf('tasks/task_%03d_%03d',j,i),'task1');
    end
%     
% end
%%
% z1=[];
% z2=[];
% figure;hold on;
% for j = 3:3
%     for i=1:length(tss)
%         if(exist(sprintf('tasks/task_%03d_%03d.mat',j,i),'file') ~= 2)
%             continue
%         end
%         load(sprintf('tasks/task_%03d_%03d.mat',j,i));
%         tshift = tss(i)*1e-3/tcoef;
%         tsind = ceil(tshift/0.2);
%         zz = gather(task1.history.N1-task1.history.N2)/NN;
%         plot(real(zz(1:end)));
%         z1(i,j) = gather(task1.history.N1(tsind)-task1.history.N2(tsind))/NN;
%         z2(i,j) = gather(sum(task1.history.N1(tsind:end)-task1.history.N2(tsind:end)))/NN/(2000-tsind+1);
%     end
% end