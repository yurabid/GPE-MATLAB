
global tcoef Ub omegaBar initphase ub1;
%   Box size
xmax = 15;
zmax = 7.5;
%   Grid size
N=128;
Nz=32;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);

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
task.Ntotal = 2000; 
task.Vtd = @TDPot;
task.gamma = 0.000;
task.user_callback = @post_process;
% task.decay_rate = 1.0e50;

% zz0 = zeros(1,20);

%%
for jj = 55:65
    ub1 = 0.025/100*jj;
    tic
    phi = task.groundstate_itp(1e-2,1e-6);
    toc
    tic
%     task.current_mu = mu(end);
%     task.current_iter = 1;
%     post_process(phi);
    task.current_iter = 0;
    task.solve_split(0.01,100,3830);
%     zz0(jj) = gather(task.history.N1(1)-task.history.N2(1))/task.Ntotal;

    toc
    save(sprintf('tasks/task_%05d',jj),'task');
end
%%
figure;
hold all;
zavg = zeros(9,200,'like',task.grid.mesh.x);
freqs = zeros(1,200,'like',task.grid.mesh.x);
angles3 = atan2(grid.mesh.x,grid.mesh.y);
curangle = 0.0;%0.125*pi;
%for jj = [90,50,18]
for jj = 58:58
    if(exist(sprintf('tasks/task_%05d.mat',jj),'file') ~= 2)
        continue
    end
    load(sprintf('tasks/task_%05d.mat',jj));
    zz=(task.history.N1-task.history.N2)./task.Ntotal;
    cl=sqrt(task.g*max(max(max(abs(task.init_state).^2.*(angles3>(-pi + curangle)).*(angles3<-curangle))))/2)/r0;
    cr=sqrt(task.g*max(max(max(abs(task.init_state).^2.*(angles3<(-pi + curangle) | angles3>-curangle))))/2)/r0;
%     zzz = sum(zz)/numel(zz);
%     fzz = abs(fft(zz-zzz));
%     freqs(jj)=(find(fzz(1:end/2)==max(fzz(round(sqrt(jj)*3):end/2))))/numel(zz);
    mu1 = grid.integrate(conj(task.init_state).*task.applyham(task.init_state).*(angles3>(-pi + curangle)).*(angles3<-curangle))/task.history.N1(1);
    mu2 = grid.integrate(conj(task.init_state).*task.applyham(task.init_state).*(angles3<(-pi + curangle) | angles3>-curangle))/task.history.N2(1);
    zavg(:,jj) = [zz(1) zz(36) zz(end) sum(zz(end-500:end))/501 task.history.dmu(36) sum(task.history.dmu(end-500:end))/501 cl cr mu1-mu2];
%     plot((1:numel(fzz))/numel(zz),fzz);
plot((1:length(zz))*tcoef,zz);
% plot((1:length(zz))*tcoef,task.history.dmu*570);
end

%%
    rhogrid=(1:256)/256*xmax;
    phigrid=(1:512)/512*2*pi+pi/2;
    drho = rhogrid(2)-rhogrid(1);
    dphi = phigrid(2)-phigrid(1);
    [RHO, PHI] = ndgrid(rhogrid,phigrid);
    YC = RHO.*sin(PHI);
    XC = RHO.*cos(PHI);
    phi_t_map = zeros(3048,512,'like',task.grid.mesh.x);
for j=1:3048
    if(exist(sprintf('snapshots/slice_%05d.mat',j),'file') ~= 2)
        continue
    end
    load(sprintf('snapshots/slice_%05d.mat',j));
    densc = interp2(grid.mesh.x2,grid.mesh.y2,densz,XC,YC,'linear');
    phi_t_map(j,:) = rhogrid*densc*drho;
end
