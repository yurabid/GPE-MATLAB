%Calculate the dynamics of GPE with the split-step method

% Initialize basic parameters
global tcoef XXg YYg;
if(exist('L','var') ~= 1)
    disp('reinitializing config parameters');
    config
end

%   Initialize grid
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
k = [ (0:N/2)*2*pi/L -(N/2-1:-1:1)*2*pi/L];
kz = [ (0:Nz/2)*2*pi/Lz -(Nz/2-1:-1:1)*2*pi/Lz];
%%
[XX,YY] = meshgrid(r,r);
XXg = gpuArray(XX);
YYg = gpuArray(YY);
if(detectCores(1) > 0 || detectCores(2) > 0)
    [XZ,ZZ] = meshgrid(r,rz);
end
[X,Y,Z] = meshgrid(r,r,rz);
[KX,KY,KZ] = meshgrid(k,k,kz);
slice=zeros(N,N);
slicez=zeros(N,Nz);

V = Vfun(X,Y,Z);
VV = V;
dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt);
clear k kz KX KY KZ;

doCoreDetect = sum(detectCores);
if(start>0)
    load('phi2.mat'); %loading the previous state
else
    load('phi.mat'); %loading the initial state
    
    if(exist('phimod','var') == 1)
        phi = phi.*phimod(X,Y,Z);
    end
    maxx = max(max(max(abs(phi))));
    if(saveSlices(1) > 0)
        slicez(:,:) = phi(N/2,:,:);
        save(sprintf('snapshots/sliceyz_%05d',start),'slicez');
    end
    if(saveSlices(2) > 0)
        slicez(:,:) = phi(:,N/2,:);
        save(sprintf('snapshots/slicexz_%05d',start),'slicez');
    end
    if(saveSlices(3) > 0)
        slice(:,:) = phi(:,:,Nz/2);
        save(sprintf('snapshots/slice_%05d',start),'slice');
    end
    
    
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
%mu_last = mu;
%%

dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
dt_outer = ddt*niter_inner;
if(start>0)
    load('params.mat')
    mu = gpuArray(MUc(start));
    %mu_last = mu;
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
        if(useTDPot > 0)
            VV = bsxfun(@plus,V,TDPot(time2));
        end
        
        %         mu_run = mu + (mu-mu_last)/niter_inner*i;
        mu_run = mu*exp(-i*ddt/tau);
        tmp = exp((mu_run - VV - g*tmp.*conj(tmp))*dt).*tmp;
    end
    
    tmp = exp((VV - mu + g*tmp.*conj(tmp))*dt/2).*tmp;
    tmp2 = abs(tmp.*conj(tmp));
    NNN = NN0*exp(-time2/tau);
    NNgpu = sum(sum(sum(tmp2)));
    NN(j) = gather(NNgpu)*h*h*hz;
    MU(j) = gather(sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (VV+g*tmp2).*tmp2)))/NNgpu);
    HH(j) = MU(j) - gather(sum(sum(sum(g*0.5*tmp2.*tmp2)))/NNgpu);
    %mu_last = mu;
    mu = mu + gpuArray((-log(NN(j)/NNN)/dt_outer)*(1+gam^2)/(2*gam));
    MUp(j) = gather(mu_run);
    MUc(j) = gather(mu);
    phi=gather(tmp);
    LL(j) = -(imag(sum(sum(sum(conj(phi).*(X.*(-circshift(phi,2)+8*circshift(phi,1)-8*circshift(phi,-1)+circshift(phi,-2)) -...
        Y.*(-circshift(phi,[0 2])+8*circshift(phi,[0 1])-8*circshift(phi,[0 -1])+circshift(phi,[0 -2]))))))))*h*hz./(12*NN(j));
    fprintf('j = %d, t = %0.3f, L = %0.6f, E = %0.6f, mu = %0.6f \n', j, time2*tcoef, LL(j), HH(j), MU(j));
    toc
    
    if(doCoreDetect > 0)
        % core detection
        coresp = [0 0 0];
        coresm = [0 0 0];
        if(detectCores(3) > 0)
            for i=1:Nz
                [coresp1, coresm1] = detect_core(phi(:,:,i),XX,YY);
                coresp = [coresp; coresp1 ones(size(coresp1,1),1).*rz(i)];
                coresm = [coresm; coresm1 ones(size(coresm1,1),1).*rz(i)];
            end
        end
        if(detectCores(1) > 0)
            for i=1:N
                [coresp1, coresm1] = detect_core(squeeze(phi(i,:,:)),XZ,ZZ);
                coresp = [coresp; coresp1(:,2) ones(size(coresp1,1),1).*r(i) coresp1(:,1)];
                coresm = [coresm; coresm1(:,2) ones(size(coresm1,1),1).*r(i) coresm1(:,1)];
            end
        end
        if(detectCores(2) > 0)
            for i=1:N
                [coresp1, coresm1] = detect_core(squeeze(phi(:,i,:)),XZ,ZZ);
                coresp = [coresp; ones(size(coresp1,1),1).*r(i) coresp1(:,2) coresp1(:,1)];
                coresm = [coresm; ones(size(coresm1,1),1).*r(i) coresm1(:,2) coresm1(:,1)];
            end
        end
        coresp = coresp(2:end,:);
        coresm = coresm(2:end,:);
        fv = isosurface(X,Y,Z,abs(phi),0.1*maxx);
        save(sprintf('snapshots/core_%05d',j),'coresp','coresm','fv');
    end
    
    if(saveSlices(1) > 0)
        slicez(:,:) = phi(N/2,:,:);
        save(sprintf('snapshots/sliceyz_%05d',j),'slicez');
    end
    if(saveSlices(2) > 0)
        slicez(:,:) = phi(:,N/2,:);
        save(sprintf('snapshots/slicexz_%05d',j),'slicez');
    end
    if(saveSlices(3) > 0)
        slice(:,:) = phi(:,:,Nz/2);
        save(sprintf('snapshots/slice_%05d',j),'slice');
    end
    save('phi2','phi');
    save('params', 'NN' ,'maxx','HH', 'LL', 'MU', 'MUc');
end
