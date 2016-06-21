h=0.3;
hz=0.2;
N=128;
Nz=32;
g = 0.0188;
start = 297;

r=h*((1:N)-N/2-0.5);
rz=hz*((1:Nz)-Nz/2-0.5);
k = [ (0:N/2)*2*pi/(r(end)-r(1)) -(N/2-1:-1:1)*2*pi/(r(end)-r(1))];
kz = [ (0:Nz/2)*2*pi/(rz(end)-rz(1)) -(Nz/2-1:-1:1)*2*pi/(rz(end)-rz(1))];
kk = zeros(N,N,Nz);
%%
[XX,YY] = meshgrid(r,r);
aXY = (atan2(XX,YY));

V=zeros(N,N,Nz);
slice=zeros(N,N);
slicez=zeros(N,Nz);
tic
rr = 10.45;
om = 0.5*4.88^2;
for ix=1:N
    for iy=1:N
        rs=(sqrt(r(ix)^2+r(iy)^2)-rr).^2;
        V(ix,iy,:) = om*rz(:).^2 + 0.5*rs;
        kk(ix,iy,:) = (k(ix)^2+k(iy)^2+kz(:).^2)/2;
    end
end
clear rr om rs k kz;
[X,Y,Z] = meshgrid(r,r,rz);
toc
f=figure;
if(start>0)
    load('phi2.mat'); %loading the previous state
else
    load('phi_128.mat'); %loading the initial state
    % save initial state
%     phi = bsxfun(@times,phi,exp(-1i*15*aXY));

%     s = 1;
%     r0 = 5;
%     xi=0.6;
%     tic
%     for ix=1:N
%         for iy=1:N
%             for iz=1:Nz
%                 theta = angle(r(ix) - r0 + 1i*r(iy));
%                 rho = sqrt((r(ix)-r0)^2+r(iy)^2);
%                 ff = tanh(sqrt((r(ix)-r0)^2+r(iy)^2)/xi)^s*exp(-1i*s*angle(r(ix) - r0 + 1i*r(iy)));
%                 phi(ix,iy,iz) = phi(ix,iy,iz)*ff;
%             end
%         end
%     end
%     toc
    
    maxx = abs(max(max(max(phi))));
%     slicez(:,:) = phi(:,N/2,:);
    slice(:,:) = phi(:,:,Nz/2);
    save(sprintf('slice_%05d',start),'slice');
%     save(sprintf('slicez_%05d',start),'slicez');
%     image(r,r,abs(slice)/maxx*64);
%     set(gca,'YDir','normal');
%     saveas(f,sprintf('sol_%05d.png',start));
end

%%
tmp = gpuArray(phi);
aXY = gpuArray(aXY);
XX = gpuArray(XX);
YY = gpuArray(YY);
gam = 0.0015; % dissipation constant gamma
ddt = 0.01;
dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
tau = 10 / 1.29366e-3; % decay constant
NN0=6.0e5; % number of particles

mu = (abs(sum(sum(sum(abs(conj(phi).*ifftn(kk.*fftn(phi))) + (V+g*phi.*conj(phi)).*phi.*conj(phi)))))*h*h*hz)./(abs(sum(sum(sum(phi.*conj(phi)))))*h*h*hz);
mu0 = mu;
mu_last = mu;
kk = gpuArray(kk);
ekk = exp(-kk*dt);
%%
niter_inner = 100; %number of internal iterations
niter_outer = 290*8; %number of external iterations
dt_outer = ddt*niter_inner;

tz = 0.5 / 1.29366e-3;
wwstart = 2.5*2*pi*1.29366e-3;
wwstop = 0.0*2*pi*1.29366e-3;

ww = wwstart;

dpi = 2*pi;
pwidth = -12.222^2;

V = gpuArray(V);
VV = V;
VVV = gpuArray.zeros(N,N);
VV2 = gpuArray.zeros(N,N);
cmax = 1300/123;
if(start>0)
    load('params.mat')
    mu_last = MUc(start);
else
    HH = zeros(niter_outer,1);
    MU = zeros(niter_outer,1);
    MUp = zeros(niter_outer,1);
    MUc = zeros(niter_outer,1);
    NN = zeros(niter_outer,1);    
    LL = zeros(niter_outer,1);
end
% main BIG cycle starts here
fftw('planner','patient');
for j=start+1:niter_outer
    tic   
    
    time=(j-1)*dt_outer;
    tmp = exp(-(VV - mu  + g*tmp.*conj(tmp))*dt/2).*tmp;

    % main SMALL cycle starts here
    for i=1:niter_inner
        tmp = ifftn(ekk.*fftn(tmp));
        
        time2=time+(i-1)*ddt;
        if(time2<tz)
            coef = time2/tz;
        elseif(time2<tz*2)
            coef = 1;
        elseif(time2<tz*3)
            coef = 1+(2*tz-time2)/tz;
        else
            ww = wwstop;
            coef = (time2-3*tz)/tz;
        end
         cs = cos(ww*time2);
         sn = sin(ww*time2);
%          VV2 = -0.146201/(1.5^2)*(XX*sn - YY*cs).^2;
         VV2 = -0.146201*(XX*sn - YY*cs).^2;
         VVV = coef*cmax*exp(VV2).*(YY*sn+XX*cs>0);

%         VV2 = (mod(ww*time2+pi-aXY,dpi)-pi);
%         VVV = coef*cmax*exp(pwidth*VV2.*VV2.*VV2.*VV2);

        VV = bsxfun(@plus,V,VVV);
        
%         mu_run = mu + (mu-mu_last)/niter_inner*i;
        mu_run = mu*exp(i*ddt/tau);
        tmp = exp((mu_run - VV - g*tmp.*conj(tmp))*dt).*tmp;
    end

    tmp2 = tmp.*conj(tmp);
    tmp = exp((VV - mu + g*tmp2)*dt/2).*tmp;
    NNN = NN0*exp(-time/tau);
    NNgpu = sum(sum(sum(tmp2)));
    NN(j) = gather(NNgpu)*h*h*hz;
    HH(j) = gather((abs(sum(sum(sum(conj(tmp).*ifftn(kk.*fftn(tmp)) + (VV+g*0.5*tmp2).*tmp2)))))/NNgpu);
    MU(j) = gather((abs(sum(sum(sum(conj(tmp).*ifftn(kk.*fftn(tmp)) + (VV+g*tmp2).*tmp2)))))/NNgpu);
    mu_last = mu;
    mu = mu +(-log(NN(j)/NNN)/dt_outer)*(1+gam^2)/(2*gam);
    MUp(j) = mu_run;
    MUc(j) = mu;
    phi=gather(tmp);
    LL(j) = (abs(sum(sum(sum(1i/(12*h)*conj(phi).*(X.*(-circshift(phi,2)+8*circshift(phi,1)-8*circshift(phi,-1)+circshift(phi,-2)) - Y.*(-circshift(phi,[0 2])+8*circshift(phi,[0 1])-8*circshift(phi,[0 -1])+circshift(phi,[0 -2]))))))))*h*h*hz./NN(j);
    fprintf('j = %d, t = %0.3f, L = %0.6f, E = %0.6f, mu = %0.6f \n', j, time2*1.29366e-3, LL(j), HH(j), MU(j));
    toc
        % core detection
%     tic 
%      coresp = [0 0 0];
%      coresm = [0 0 0];
% 
%     for i=1:Nz
%         slice(:,:) = phi(:,:,i);
%         cores = detect_core(slice,h);
%         [y,x] = find(cores>0);
%         ncores = length(x);
%         for ii = 1:ncores
%             coresp = [coresp; r(x(ii))-h/2 r(y(ii))-h/2 rz(i)];
%         end
%         [y,x] = find(cores<0);
%         ncores = length(x);
%         for ii = 1:ncores
%             coresm = [coresm; r(x(ii))-h/2 r(y(ii))-h/2 rz(i)];
%         end
%     end    
%     coresp = coresp(2:end,:);
%     coresm = coresm(2:end,:);
%     toc
%      fv = isosurface(X,Y,Z,abs(phi),0.1*maxx);
%      fv2 = isosurface(X,Y,Z,abs(phi),0.4*maxx);    
%      save(sprintf('core_%05d',j),'coresp','coresm','fv','fv2','VVV');
    
%     slicez(:,:) = phi(:,N/2,:);
    slice(:,:) = phi(:,:,Nz/2);
    
% [sx,sy] = gradient(slice);
% [sxc,syc] = gradient(conj(slice));
% vx = -imag(slice.*sxc - conj(slice).*sx)./(abs(slice).^2).*(abs(slice)>1e-4)./(2*h);
% vy = -imag(slice.*syc - conj(slice).*sy)./(abs(slice).^2).*(abs(slice)>1e-4)./(2*h);
    
%     image(r,r,abs(slice)/maxx*64);
%     image(r,r,sqrt(vx.^2+vy.^2)./sqrt(mu0/2)*64);
%     set(gca,'YDir','normal');
%     saveas(f,sprintf('sol_%05d.png',j));
%     save(sprintf('slicez_%05d',j),'slicez');
    save(sprintf('slice_%05d',j),'slice');
    save('phi2','phi');
    save('params', 'NN' ,'mu','maxx','HH', 'LL', 'MU', 'MUp', 'MUc');
end
close all;