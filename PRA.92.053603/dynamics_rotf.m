% Initialize basic parameters
%   Box size
L = 26;
Lz = 26;
%   Grid size
N=128;
Nz=128;
%   Physical parameters
g = 0.0276; %0.092853;
NN0=1.0e6; % number of particles
start = 3030;
tcoef = 1/(100*2*pi); %time scale

gam = 0.03; % dissipation constant gamma
ddt = 0.007;
dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
tau = 1000000 / tcoef; % decay constant

vm = 34;
vz = 18;
alpha = 0.8;
beta = 0.2;
m = 1;

%   Initialize grid
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
k = [ (0:N/2)*2*pi/L -(N/2-1:-1:1)*2*pi/L];
kz = [ (0:Nz/2)*2*pi/Lz -(Nz/2-1:-1:1)*2*pi/Lz];
%%
[XX,YY] = meshgrid(r,r);
[X,Y,Z] = meshgrid(r,r,rz);
RS = X.^2 + Y.^2;
[KX,KY,KZ] = meshgrid(k,k,kz);
slice=zeros(N,N);
slicez=zeros(N,Nz);

%V = gpuArray(0.5*Z.^2 + 0.5*X.^2+0.5*Y.^2);
V = gpuArray(0.5*Z.^2 * 9 + 0.5*RS + vm*beta^(2*m)*RS.^m.*exp(-m*(beta^2*RS - 1)) + vz*exp(-alpha^2*Z.^2));
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt);
clear k kz KX KY KZ;

%f=figure;
if(start>0)
    load('phi2.mat'); %loading the previous state
else
    load('phi.mat'); %loading the initial state
	
    s = 1;
    sr = 1;
    r0 = 0.5;
    rr0 = 6;
    z0 = 0.1;
    xi=0.2;

    phi = phi.*tanh(sqrt((X-r0).^2+Y.^2)/xi).^s.*exp(-1i*s*atan2(X - r0,Y)); % imprint s-charged off-center vortex
    phi = phi.*tanh(sqrt((sqrt(RS)-rr0).^2+(Z-z0).^2)/xi).^sr.*exp(-1i*sr*atan2(sqrt(RS)-rr0,Z-z0)); % imprint vortex ring
    maxx = max(max(max(abs(phi))));
    slice(:,:) = phi(:,:,Nz/2);
    save(sprintf('snapshots/slice_%05d',start),'slice');
    slicez(:,:) = phi(:,N/2,:);
    save(sprintf('snapshots/slicez_%05d',start),'slicez');
end
clear RS;

%%
tmp = gpuArray(phi);
Xg = gpuArray(X);
Yg = gpuArray(Y);
%aXY = atan2(XXg,YYg);
%clear X Y Z;
tmp2 = abs(tmp.*conj(tmp));
mu = (abs(sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (V+g*tmp2).*tmp2))))*h*h*hz)./(abs(sum(sum(sum(tmp2))))*h*h*hz);
mu_last = mu;
%%
niter_inner = 100; %number of internal iterations
niter_outer = 5000; %number of external iterations
dt_outer = ddt*niter_inner;

omega = 0.25;
n_cn = 2;


%VV = V;
VVV = gpuArray.zeros(N,N);

if(start>0)
    load('params.mat')
    mu_last = gpuArray(MUc(start));
    mu = mu_last;
else
    HH = zeros(niter_outer,1);
    MU = zeros(niter_outer,1);
    MUp = zeros(niter_outer,1);
    MUc = zeros(niter_outer,1);
    NN = zeros(niter_outer,1);    
    LL = zeros(niter_outer,1);
end
% main BIG cycle starts here
for j=start+1:niter_outer
    tic   
    
    time=(j-1)*dt_outer;
    tmp = exp(-(V - mu  + g*tmp2)*dt/2).*tmp;

    % main SMALL cycle starts here
    for i=1:niter_inner
        tmp = ifftn(ekk.*fftn(tmp));
        
        if(omega ~= 0)
            lphi = tmp;
            for ii = 1:n_cn
                lphi = tmp + dt*1i*omega/(12*h)*(Xg.*(-circshift(lphi,[2 0])+8*circshift(lphi,[1 0])-8*circshift(lphi,[-1 0])+circshift(lphi,[-2 0])) - ...
                    Yg.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
                lphi = 0.5*(tmp+lphi);
            end
            tmp = tmp + dt*1i*omega/(12*h)*(Xg.*(-circshift(lphi,2)+8*circshift(lphi,1)-8*circshift(lphi,-1)+circshift(lphi,-2)) - ...
                Yg.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
        end
		
		
		
        time2=time+(i-1)*ddt;
        mu_run = mu*exp((i-1)*ddt/tau);
        tmp = exp((mu_run - V - g*tmp.*conj(tmp))*dt).*tmp;
    end

    tmp = exp((V - mu + g*tmp.*conj(tmp))*dt/2).*tmp;
	tmp2 = abs(tmp.*conj(tmp));
    NNN = NN0*exp(-time2/tau);
    NNgpu = sum(sum(sum(tmp2)));
    NN(j) = gather(NNgpu)*h*h*hz;
    MU(j) = gather(sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (V+g*tmp2).*tmp2)))/NNgpu);
    HH(j) = MU(j) - gather(sum(sum(sum(g*0.5*tmp2.*tmp2)))/NNgpu);
    LL(j) = gather((abs(sum(sum(sum(conj(tmp).*(Xg.*(-circshift(tmp,[2 0])+8*circshift(tmp,[1 0])-8*circshift(tmp,[-1 0])+circshift(tmp,[-2 0])) -...
												   Yg.*(-circshift(tmp,[0 2])+8*circshift(tmp,[0 1])-8*circshift(tmp,[0 -1])+circshift(tmp,[0 -2]))))))))./(12*NNgpu*h));
    mu_last = mu;
    mu = mu + gpuArray((-log(NN(j)/NNN)/dt_outer)*(1+gam^2)/(2*gam));
    MUp(j) = gather(mu_run);
    MUc(j) = gather(mu);
    phi=gather(tmp);
    fprintf('j = %d, t = %0.3f, L = %0.6f, E = %0.6f, mu = %0.6f \n', j, time2*tcoef, LL(j), HH(j), MU(j));
    toc
        
    % core detection
    tic 
    coresp = [0 0 0];
    coresm = [0 0 0];
    for i=1:Nz
        [coresp1, coresm1] = detect_core(phi(:,:,i),XX,YY);
        coresp = [coresp; coresp1 ones(size(coresp1,1),1).*rz(i)];
        coresm = [coresm; coresm1 ones(size(coresm1,1),1).*rz(i)];
    end   
    for i=1:N
        [coresp1, coresm1] = detect_core(squeeze(phi(i,:,:)),XX,YY);
        coresp = [coresp; coresp1(:,2) ones(size(coresp1,1),1).*r(i) coresp1(:,1)];
        coresm = [coresm; coresm1(:,2) ones(size(coresm1,1),1).*r(i) coresm1(:,1)];
    end 
     for i=1:N
         [coresp1, coresm1] = detect_core(squeeze(phi(:,i,:)),XX,YY);
         coresp = [coresp; ones(size(coresp1,1),1).*r(i) coresp1(:,2) coresp1(:,1)];
         coresm = [coresm; ones(size(coresm1,1),1).*r(i) coresm1(:,2) coresm1(:,1)];
     end 
    coresp = coresp(2:end,:);
    coresm = coresm(2:end,:);
    toc
    fv = isosurface(X,Y,Z,abs(phi),0.1*maxx);    
    save(sprintf('snapshots/core_%05d',j),'coresp','coresm','fv');
    
    slice(:,:) = phi(:,:,Nz/2);
    save(sprintf('snapshots/slice_%05d',j),'slice');
    slicez(:,:) = phi(:,N/2,:);
    save(sprintf('snapshots/slicez_%05d',j),'slicez');
    save('phi2','phi');
    save('params', 'NN' ,'maxx','HH', 'LL', 'MU', 'MUp', 'MUc');
end
close all;