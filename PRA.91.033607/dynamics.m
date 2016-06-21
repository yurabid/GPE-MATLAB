h=0.15;
hz=0.2;
N=256;
Nz=32;
g = 0.0188;
start = 1854;

gam = 0.0015; % dissipation constant gamma
ddt = 0.005;
dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
tau = 10000 / 1.29366e-3; % decay constant
NN0=6.0e5; % number of particles

r=h*((1:N)-N/2-0.5);
rz=hz*((1:Nz)-Nz/2-0.5);
k = [ (0:N/2)*2*pi/(r(end)-r(1)) -(N/2-1:-1:1)*2*pi/(r(end)-r(1))];
kz = [ (0:Nz/2)*2*pi/(rz(end)-rz(1)) -(Nz/2-1:-1:1)*2*pi/(rz(end)-rz(1))];
%%
[XX,YY] = meshgrid(r,r);
[X,Y,Z] = meshgrid(r,r,rz);
[KX,KY,KZ] = meshgrid(k,k,kz);
slice=zeros(N,N);
slicez=zeros(N,Nz);

r0 = 10.45;
om = 0.5*4.88^2;
V = gpuArray(om*Z.^2 + 0.5*(sqrt(X.^2+Y.^2)-r0).^2);
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt);
clear k kz KX KY KZ;

%f=figure;
if(start>0)
    load('phi2.mat'); %loading the previous state
else
    load('phi256_1.mat'); %loading the initial state
	
%     phi = bsxfun(@times,phi,exp(-1i*15*atan2(XX,YY))); % imprint 15-charged vortex in the enter
%     s = 1;
%     r0 = 5;
%     xi=0.6;
%     phi = bsxfun(@times,phi,tanh(sqrt((XX-r0).^2+YY^2)/xi).^s*exp(-1i*s*atan2(XX - r0,YY))); % imprint s-charged off-center vortex
    
    maxx = max(max(max(abs(phi))));
    slice(:,:) = phi(:,:,Nz/2);
    save(sprintf('snapshots/slice_%05d',start),'slice');
end

%%
tmp = gpuArray(phi);
XXg = gpuArray(XX);
YYg = gpuArray(YY);
%aXY = atan2(XXg,YYg);

tmp2 = abs(tmp.*conj(tmp));
mu = (abs(sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (V+g*tmp2).*tmp2))))*h*h*hz)./(abs(sum(sum(sum(tmp2))))*h*h*hz);
mu_last = mu;
%%
niter_inner = 40; %number of internal iterations
niter_outer = 5000; %number of external iterations
dt_outer = ddt*niter_inner;

tz = 1 / 1.29366e-3;
wwstart = -1.5*2*pi*1.29366e-3;
wwstop = 0.0*2*pi*1.29366e-3;
ww = wwstart;

VV = V;
VVV = gpuArray.zeros(N,N);
cmax = 1300/123;
if(start>0)
    load('params.mat')
    mu_last = gpuArray(MUc(start));
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
    tmp = exp(-(VV - mu  + g*tmp2)*dt/2).*tmp;

    % main SMALL cycle starts here
    for i=1:niter_inner
        tmp = ifftn(ekk.*fftn(tmp));
        
        time2=time+(i-1)*ddt;
%         if(time2<tz)
            coef = time2/tz;
%         elseif(time2<tz*2)
%             coef = 1;
%         elseif(time2<tz*3)
%             coef = 1+(2*tz-time2)/tz;
%         else
%             %ww = wwstop;
%             %coef = (time2-3*tz)/tz;
%             coef = 0;
%         end
        cs = cos(ww*time2);
        sn = sin(ww*time2);
        VVV = coef*cmax*exp(-0.146201*(XXg*sn - YYg*cs).^2).*(YYg*sn+XXg*cs>0);
        VV = bsxfun(@plus,V,VVV);
        
%         mu_run = mu + (mu-mu_last)/niter_inner*i;
        mu_run = mu*exp(i*ddt/tau);
        tmp = exp((mu_run - VV - g*tmp.*conj(tmp))*dt).*tmp;
    end

    tmp = exp((VV - mu + g*tmp.*conj(tmp))*dt/2).*tmp;
	tmp2 = abs(tmp.*conj(tmp));
    NNN = NN0*exp(-time2/tau);
    NNgpu = sum(sum(sum(tmp2)));
    NN(j) = gather(NNgpu)*h*h*hz;
    MU(j) = gather(sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (VV+g*tmp2).*tmp2)))/NNgpu);
    HH(j) = MU(j) - gather(sum(sum(sum(g*0.5*tmp2.*tmp2)))/NNgpu);
    mu_last = mu;
    mu = mu + gpuArray((-log(NN(j)/NNN)/dt_outer)*(1+gam^2)/(2*gam));
    MUp(j) = gather(mu_run);
    MUc(j) = gather(mu);
    phi=gather(tmp);
    LL(j) = (abs(sum(sum(sum(conj(phi).*(X.*(-circshift(phi,2)+8*circshift(phi,1)-8*circshift(phi,-1)+circshift(phi,-2)) -...
												   Y.*(-circshift(phi,[0 2])+8*circshift(phi,[0 1])-8*circshift(phi,[0 -1])+circshift(phi,[0 -2]))))))))*h*hz./(12*NN(j));
    fprintf('j = %d, t = %0.3f, L = %0.6f, E = %0.6f, mu = %0.6f \n', j, time2*1.29366e-3, LL(j), HH(j), MU(j));
    toc
        % core detection
%    tic 
%    coresp = [0 0 0];
%    coresm = [0 0 0];
%    for i=1:Nz
%        [coresp1, coresm1] = detect_core(phi(:,:,i),XX,YY);
%        coresp = [coresp; coresp1 ones(size(coresp1,1),1).*rz(i)];
%        coresm = [coresm; coresm1 ones(size(coresm1,1),1).*rz(i)];
%    end    
%    coresp = coresp(2:end,:)
%    coresm = coresm(2:end,:);
%    toc
%    fv = isosurface(X,Y,Z,abs(phi),0.1*maxx);
%    fv2 = isosurface(X,Y,Z,abs(phi),0.4*maxx);    
%    save(sprintf('core_%05d',j),'coresp','coresm','fv');
    
    slice(:,:) = phi(:,:,Nz/2);
    
   
%     image(r,r,abs(slice)/maxx*64);
%     set(gca,'YDir','normal');
%     saveas(f,sprintf('sol_%05d.png',j));
    save(sprintf('snapshots/slice_%05d',j),'slice');
    save('phi2','phi');
    save('params', 'NN' ,'maxx','HH', 'LL', 'MU', 'MUp', 'MUc');
end
close all;