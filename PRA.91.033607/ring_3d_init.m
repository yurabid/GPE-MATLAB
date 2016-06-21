h=0.15;
hz=0.2;
N=256;
Nz=32;
g = 0.0188;
r=h*((1:N)-N/2-0.5);
rz=hz*((1:Nz)-Nz/2-0.5);
k = [ (0:N/2)*2*pi/(r(end)-r(1)) -(N/2-1:-1:1)*2*pi/(r(end)-r(1))];
kz = [ (0:Nz/2)*2*pi/(rz(end)-rz(1)) -(Nz/2-1:-1:1)*2*pi/(rz(end)-rz(1))];
kk = zeros(N,N,Nz);
V=zeros(N,N,Nz);
phi0=zeros(N,N,Nz);
[X,Y,Z] = meshgrid(r,r,rz);

%% Construct the potential V, initial condition phi0 and momentum array kk
tic
rr = 10.45;
om = 0.5*4.88^2;

for ix=1:N
    for iy=1:N
        rs=r(ix)^2+r(iy)^2;
        V(ix,iy,:) = om*rz(:).^2 + 0.5*(sqrt(rs)-rr).^2;
        phi0(ix,iy,:) = exp(-rs/(20)-rz(:).^2/(20));
        kk(ix,iy,:) = (k(ix).^2+k(iy).^2+kz(:).^2)/2;
    end
end
toc

%%
dt = 0.004; % time step
ekk = exp(-kk*dt);

NN=sqrt(6.0e5); % square root of number of particles
niter = 500; % number of ITP itrerations
phi = phi0*NN/(sqrt(sum(sum(sum(abs(phi0).^2)))*h*h*hz));
tmp = gpuArray(phi);
V = gpuArray(V);
MU = zeros(niter,1);
MU2 = zeros(niter,1);
clear phi0 Vz;
fftw('planner', 'patient');
tic
for i=1:niter
    tmp = exp(-(V + g*tmp.*conj(tmp))*dt*0.5).*tmp;
    tmp = ifftn(ekk.*fftn(tmp));
    tmp = exp(-(V + g*tmp.*conj(tmp))*dt*0.5).*tmp;
    mu = NN/sqrt(sum(sum(sum(tmp.*conj(tmp))))*h*h*hz);
    tmp=tmp*mu;
    MU(i) = gather(1/dt * log(mu));
    MU2(i) = gather((sum(sum(sum(abs(conj(tmp).*ifftn(kk.*fftn(tmp))) + (V+g*tmp.*conj(tmp)).*tmp.*conj(tmp) ))))*h*h*hz./NN^2);
end
phi=gather(tmp);
toc
%% Save the results
save('MU','MU','MU2');
save('phi','phi');