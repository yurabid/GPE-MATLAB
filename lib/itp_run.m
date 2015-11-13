%Calculate the stationary state of GPE with Imaginary Time Propagation method
%
% Initialize basic parameters
global XXg YYg;
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
[XX,YY] = meshgrid(r,r);
XXg = gpuArray(XX);
YYg = gpuArray(YY);
if(useTDPot > 0)
    V = bsxfun(@plus,V,TDPot(0));
end
if(exist('phi0','var') ~= 1)
    disp('no phi0 found, generating default');
    phi0 = exp((X.^2+Y.^2+Z.^2)/min(L,Lz)^2);
end
phi = phi0*sqrt(NN0/(sum(sum(sum(abs(phi0).^2)))*h*h*hz));

kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt_itp);
clear phi0 Z KX KY KZ XXg YYg XX YY;
MU = zeros(niter,1,'gpuArray');
MU2 = zeros(niter,1,'gpuArray');
tmp = gpuArray(phi);

%% Do the main calculation
tic
tmp2 = abs(tmp.*conj(tmp));
for i=1:niter
    tmp = exp(-(V + g*tmp2)*dt_itp*0.5).*tmp;
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
    tmp = exp(-(V + g*tmp.*conj(tmp))*dt_itp*0.5).*tmp;
	tmp2 = abs(tmp.*conj(tmp));
    mu = sqrt(NN0/(sum(sum(sum(tmp2)))*h*h*hz));
    tmp=tmp*mu;
	tmp2 = tmp2*mu^2;
    MU(i) = mu;
    MU2(i) = sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (V+g*tmp2).*tmp2)));
end
toc

%% Gather results and save
MU = 1/dt_itp * log(gather(MU));
MU2 = gather(MU2)*h*h*hz./NN0;
phi=gather(tmp);
save('MU','MU','MU2');
save('phi','phi');
%imagesc(r,r,abs(phi(:,:,Nz/2)));