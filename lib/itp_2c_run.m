%Calculate the stationary state of GPE with Imaginary Time Propagation method
%
% Initialize basic parameters

if(exist('L','var') ~= 1)
    disp('reinitializing config parameters');
    config
end
Nstates = 2;
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
cphi = cell(1,Nstates);
cphi2 = cell(1,Nstates);
M = gpuArray(zeros(Nstates,Nstates));
rng('shuffle');
load('phi.mat');
for i=1:Nstates
%     phi0 = rand(N,N,Nz)+1i*rand(N,N,Nz);
    phi0 = phi;
    cphi{i} = gpuArray(phi0*sqrt(1/(sum(sum(sum(abs(phi0).^2)))*h*h*hz)));
end
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt_itp);
clear phi0 Z KX KY KZ;
MU = zeros(niter,Nstates,'gpuArray');
MU2 = zeros(niter,1,'gpuArray');
% tmp = gpuArray(phi);

%% Do the main calculation
tic
% tmp2 = abs(tmp.*conj(tmp));
for i=1:niter
    tmp1 = cphi{1};
    %for j=1:2
        cphi{1} = exp(-(V + g*NN0*(0.75*cphi{1}.*conj(cphi{1}) + 0.5*cphi{2}.*conj(cphi{2})))*dt_itp*0.5).*cphi{1};
        cphi{1} = ifftn(ekk.*fftn(cphi{1}));
        if(omega ~= 0)
            lphi = cphi{1};
            for ii = 1:n_cn
                lphi = cphi{1} - 0.25*g*NN0*dt_itp*conj(lphi).*cphi{2}.^2 + dt_itp*1i*omega/(12*h)*(X.*(-circshift(lphi,[2 0])+8*circshift(lphi,[1 0])-8*circshift(lphi,[-1 0])+circshift(lphi,[-2 0])) - ...
                    Y.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
                lphi = 0.5*(cphi{1}+lphi);
            end
            cphi{1} = cphi{1} - 0.25*g*NN0*dt_itp*conj(lphi).*cphi{2}.^2 + dt_itp*1i*omega/(12*h)*(X.*(-circshift(lphi,2)+8*circshift(lphi,1)-8*circshift(lphi,-1)+circshift(lphi,-2)) - ...
                Y.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
        end
        cphi{1} = exp(-(V + g*NN0*(0.75*cphi{1}.*conj(cphi{1}) + 0.5*cphi{2}.*conj(cphi{2})))*dt_itp*0.5).*cphi{1};

        cphi{2} = exp(-(V + g*NN0*(0.5*tmp1.*conj(tmp1) + 0.75*cphi{2}.*conj(cphi{2})))*dt_itp*0.5).*cphi{2};
        cphi{2} = ifftn(ekk.*fftn(cphi{2}));
        if(omega ~= 0)
            lphi = cphi{2};
            for ii = 1:n_cn
                lphi = cphi{2} - 0.25*g*NN0*dt_itp*conj(lphi).*tmp1.^2 + dt_itp*1i*omega/(12*h)*(X.*(-circshift(lphi,[2 0])+8*circshift(lphi,[1 0])-8*circshift(lphi,[-1 0])+circshift(lphi,[-2 0])) - ...
                    Y.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
                lphi = 0.5*(cphi{2}+lphi);
            end
            cphi{2} = cphi{2} - 0.25*g*NN0*dt_itp*conj(lphi).*tmp1.^2 + dt_itp*1i*omega/(12*h)*(X.*(-circshift(lphi,2)+8*circshift(lphi,1)-8*circshift(lphi,-1)+circshift(lphi,-2)) - ...
                Y.*(-circshift(lphi,[0 2])+8*circshift(lphi,[0 1])-8*circshift(lphi,[0 -1])+circshift(lphi,[0 -2])));
        end
        cphi{2} = exp(-(V + g*NN0*(0.5*tmp1.*conj(tmp1) + 0.75*cphi{2}.*conj(cphi{2})))*dt_itp*0.5).*cphi{2};
    %end
    for j=1:Nstates
%         for jj=1:j
%             M(jj,j) = sum(sum(sum(cphi{j}.*conj(cphi{jj}))));
%             M(j,jj) = conj(M(jj,j));
%         end
        mu = sum(sum(sum(cphi{j}.*conj(cphi{j}))))*h*h*hz;
        MU(i,j) = mu;
        cphi{j} = cphi{j}/sqrt(mu);
    end
%      relphase = angle(sum(sum(sum(cphi{1}.*conj(cphi{2})))));
%      cphi{2} = cphi{2} * exp(1i*relphase);
    %[U,S] = eig(M);
%     mu = diag(M)*h*h*hz/1;
%     for j=1:Nstates
        %cphi2{j} = gpuArray(zeros(N,N,Nz));
        %for jj=1:Nstates
%             cphi{j} = 1/sqrt(mu(j))*cphi{jj};
        %end
%     end
    %cphi = cphi2;
% 	tmp2 = abs(tmp.*conj(tmp));
%     mu = sqrt(NN0/(sum(sum(sum(tmp2)))*h*h*hz));
%     tmp=tmp*mu;
% 	tmp2 = tmp2*mu^2;
%     MU(i,:) = mu.';
%     MU2(i) = M(1,2);
%    MU2(i) = sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp)) + (V+g*tmp2).*tmp2))));
end
toc

%% Gather results and save
MU = -1/dt_itp * log(gather(MU));
% MU2 = gather(MU2)*h*h*hz./NN0;
for i=1:Nstates
    cphi2{i}=gather(cphi{i});
end
save('MU','MU','MU2');
save('cphi','cphi2');
% imagesc(r,r,abs(phi(:,:,Nz/2)));
