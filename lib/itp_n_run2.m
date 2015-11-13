%Calculate the stationary state of GPE with Imaginary Time Propagation method
%
% Initialize basic parameters

if(exist('L','var') ~= 1)
    disp('reinitializing config parameters');
    config
end
Nstates = 4;

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

for i=1:Nstates
    phi0 = rand(N,N,Nz)+1i*rand(N,N,Nz);
    cphi{i} = gpuArray(phi0*sqrt(NN0/(sum(sum(sum(abs(phi0).^2)))*h*h*hz)));
end
kk = gpuArray((KX.^2+KY.^2+KZ.^2)/2);
ekk = exp(-kk*dt_itp);
clear phi0 X Y Z KX KY KZ;
MU = zeros(100,Nstates,'gpuArray');
MU2 = zeros(100,Nstates,'gpuArray');
% tmp = gpuArray(phi);

%% Do the main calculation

% tmp2 = abs(tmp.*conj(tmp));
for jjj=1:Nstates
    tic
%     M = gpuArray(zeros(jjj,jjj));
    differ = 1;
    differ2 = 1;
    mu=1;
    i = 1;
    while differ > 1.0e-9 && differ2 > 1.0e-9
        cphi{jjj} = exp(-(V + g*cphi{jjj}.*conj(cphi{jjj}))*dt_itp*0.5).*cphi{jjj};
        cphi{jjj} = ifftn(ekk.*fftn(cphi{jjj}));
        cphi{jjj} = exp(-(V + g*cphi{jjj}.*conj(cphi{jjj}))*dt_itp*0.5).*cphi{jjj};
        if(i<300)
            cphi{jjj} = cphi{jjj} + (rand(N,N,Nz)+1i*rand(N,N,Nz));
        end
        for j=1:jjj-1
%             for jj=1:j
%                 M(jj,j) = sum(sum(sum(cphi{j}.*conj(cphi{jj}))));
%                 M(j,jj) = conj(M(jj,j));
%             end
            cphi{jjj} = cphi{jjj} - cphi{j}*sum(sum(sum(cphi{jjj}.*conj(cphi{j}))))*h*h*hz/NN0;
        end

%         [U,S] = eig(M);
%         mu = diag(S)*h*h*hz/NN0;
% %         for j=1:Nstates
%             tmp = gpuArray(zeros(N,N,Nz));
%             for jj=1:jjj
%                 tmp = tmp + U(jj,j)/sqrt(mu(j))*cphi{jj};
%             end
% %         end
%         cphi{jjj} = tmp;
     	tmp2 = abs(cphi{jjj}.*conj(cphi{jjj}));
        mu = sqrt(NN0/(sum(sum(sum(tmp2)))*h*h*hz));
        cphi{jjj}=cphi{jjj}*mu;
     	tmp2 = tmp2*mu^2;
        MU(i,jjj) = mu;
        MU2(i,jjj) = sum(sum(sum(abs(conj(cphi{jjj}).*ifftn(kk.*fftn(cphi{jjj})) + (V+g*tmp2).*tmp2))));
        if(i>1) 
            differ = abs((MU(i,jjj)-MU(i-1,jjj))/MU(i,jjj));
        end
        if(i>100) 
            differ2 = abs((MU(i,jjj)-MU(i-100,jjj))/MU(i,jjj));
        end
        i = i + 1;
    end
    i
    toc
end


%% Gather results and save
MU = 1/dt_itp * log(gather(MU));
% MU2 = gather(MU2)*h*h*hz./NN0;
for i=1:Nstates
    cphi2{i}=gather(cphi{i});
end
save('MU','MU','MU2');
save('phi','cphi2');
% imagesc(r,r,abs(phi(:,:,Nz/2)));