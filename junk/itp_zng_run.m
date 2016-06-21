%Calculate the stationary state of GPE with Imaginary Time Propagation method
%
% Initialize basic parameters

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

%% Initialize the potential V, initial condition phi0, momentum array kk and other necessary arrays
[X,Y,Z] = meshgrid(r,r,rz);
[KX,KY,KZ] = meshgrid(k,k,kz);
V = Vfun(X,Y,Z);
% phi0 = exp((X.^2+Y.^2+Z.^2)/min(L,Lz)^2);
phi0 = rand(N,N,Nz);
phi = phi0*sqrt(NN0/(sum(sum(sum(abs(phi0).^2)))*h*h*hz));
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt_itp);
clear phi0 Z KX KY KZ;
if(omega == 0)
    clear X Y;
end
MU = zeros(niter,1,'gpuArray');
KE = zeros(niter,1,'gpuArray');
PE = zeros(niter,1,'gpuArray');
tmp = gpuArray(phi);
nt = tmp*0;
NNN = NN0;
NNt = 0;

%% Do the main calculation
tic
tmp2 = abs(tmp.*conj(tmp));
for i=1:niter
    tmp = exp(-(V + g*(tmp2+2*nt))*dt_itp*0.5).*tmp;
    tmp = ifftn(ekk.*fftn(tmp));
    if(omega ~= 0)
        lphi = tmp;
        for ii = 1:n_cn
            lphi = tmp + dt_itp*1i*omega/(12*h)*(X.*(-circshift(lphi,[2 0])+8*circshift(lphi,[1 0])-8*circshift(lphi,[-1 0])+circshift(lphi,[-2 0])) - ...
                Y.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
            lphi = 0.5*(tmp+lphi);
        end
        tmp = tmp + dt_itp*1i*omega/(12*h)*(X.*(-circshift(lphi,2)+8*circshift(lphi,1)-8*circshift(lphi,-1)+circshift(lphi,-2)) - ...
            Y.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
    end
    tmp = exp(-(V + g*(tmp.*conj(tmp)+2*nt))*dt_itp*0.5).*tmp;
    tmp2 = abs(tmp.*conj(tmp));
    NN = sum(sum(sum(tmp2)));
    if(TT>0 && i>niter/2) % for better performance and stability we do one half of total iterations without a thermal cloud
        NNt = sum(sum(sum(nt)))*h*h*hz;
        NNN = NN0 - NNt;
        if(NNN<1)
            NNN=1;
            nt = nt*NN0/NNt; % we need to get the correct total number of particles even above Tc
        end
    end
    mu = sqrt(NNN/(NN*h*h*hz));
    tmp=tmp*mu;
    tmp2 = tmp2*mu^2;
    MU(i) = mu;
    KE(i) = sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))))))/(NN*mu^2);
    PE(i) = sum(sum(sum((V+g*(tmp2+2*nt)).*tmp2)))/(NN*mu^2);
%     MU2(i) = KE(i)+PE(i);
    if(TT>0 && i>niter/2)
        mmu = min(MU2(i),min(min(min(V+2*g*(tmp2+nt))))-1e-10); % compensate for possibly inaccurate chem. pot. calculation
        ntt=((TT/2/pi)^(3/2)*polylog(3/2,exp((mmu-V-2*g*(tmp2+nt))/TT))); % averaging increases stability for high temperatures
        if(NNN<=1)
            NNt = sum(sum(sum(ntt)))*h*h*hz;
            nt = (nt+ntt*NN0/NNt)*0.5;
        else
            nt=(nt+ntt)*0.5;
        end
    end
end
toc

%% Gather results and save
MU = 1/dt_itp * log(gather(MU));
KE = gather(KE);
PE = gather(PE);
MU2 = KE+PE;
phi=gather(tmp);
save('MU','MU','MU2');
save('phi','phi');
if(TT>0)
    nt=gather(nt);
    save('therm','nt');
end
%imagesc(r,r,abs(phi(:,:,Nz/2)));