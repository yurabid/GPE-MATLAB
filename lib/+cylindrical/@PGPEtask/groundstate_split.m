function [phi, varargout] = groundstate_split(task,dt,eps,phi0)
% groundstate_split - Calculate the stationary state of PGPE with split step Imaginary Time Propagation method.
%
%  Usage :
%    phi = task.groundstate_itp(dt,eps)
%    phi = task.groundstate_itp(dt,eps,phi0)
%    [phi, mu] = task.groundstate_itp(dt,eps)
%    [phi, mu, mu2] = task.groundstate_itp(dt,eps)
%  Input
%    dt    :  evolution time step
%    eps   :  desired accuracy (applied to chemical potential)
%    phi0  :  initial approximation of the wave function
%  Output
%    phi      :  calculated stationary state
%    mu       :  array of chemical potential values from norm decrease
%    mu2      :  array of chemical potential from integral evaluation

grid = task.grid;
V = task.getVtotal(0);
task.V0 = V;
g = task.g;%*task.Ntotal;
if(nargin <= 3)
    phi = grid.normalize(rand(size(grid.etot),'like',V) + 1i*rand(size(grid.etot),'like',V));
else
    phi = grid.normalize(phi0);
end
% etot = grid.etot+(grid.etotphi-1/8)/grid.r0^2; 
ekk = exp(-grid.etot*dt*0.5);
ekkphi = exp(-grid.rdphimat*dt*0.5);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);
i = 0;
iswitch = 50;
phirzp = grid.sp2grid(phi);
tmp = abs(phirzp).^2;
while true
    
    i=i+1;
    phi = ekk.*phi;
    phir = grid.ifftr(phi);
    phir = ekkphi.*phir;
    phirzp = grid.ifftzphi(phir);
    phirzp = exp(-(V+g*tmp)*dt).*phirzp;
    phir = grid.fftzphi(phirzp);
    phir = ekkphi.*phir;
    phi = grid.fftr(phir);
%     phi = grid.grid2sp(phirzp);
    phi = ekk.*phi;
    phirzp = grid.sp2grid(phi);
%     phir = grid.ifftr(phi);
%     phir = ekkphi.*phir;
%     phirzp = grid.ifftzphi(phir);
%     phirzp = exp(-(V+g*tmp)*dt*0.5).*phirzp;
%     phi = grid.grid2sp(phirzp);

    if(task.mu_init > 0)
        mu = exp(task.mu_init*dt);
        MU(i) = sum(abs(grid.to1d(phi)).^2)*mu^2;
    else
        mu = sqrt(task.Ntotal/sum(abs(grid.to1d(phi)).^2));
        MU(i) = log(mu)/dt;
    end    
    
    phi=phi*mu;
    ncur = sum(abs(grid.to1d(phi)).^2);
    phirzp = phirzp*mu;%grid.sp2grid(phi);
    tmp = abs(phirzp).^2;
    
    if(nargout >= 3)
        h0 = real(sum(grid.to1d(conj(phi).*(grid.applyh0(phi) + grid.grid2sp(V.*phirzp)))));
        h1 = grid.integrate_grid(g*abs(phirzp).^4);
        MU2(i) = (h0 + h1)/ncur;
        EE(i) = (h0 + h1/2)/ncur;
    end
    
    if((i-iswitch)>10 && mod(i,10) == 0)
%         delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt/MU(i);
        delta = abs((MU(i)-MU(i-10))^2/(MU(i)-2*MU(i-10)+MU(i-20)))/MU(i);
        if(delta < eps)
            if (dt<eps || dt<1e-5)
                break;
            else
                dt = dt/1.5;
                ekk = exp(-grid.etot*dt*0.5);
                ekkphi = exp(-grid.rdphimat*dt*0.5);
                iswitch = i;
            end
        end
    end
    
    if(i>=10000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:nnz(MU));
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = MU2(1:nnz(MU2));
    varargout{2} = MU2;
end
if(nargout >= 4)
    EE = EE(1:nnz(EE));
    varargout{3} = EE;
end
task.init_state = phi;
phi=grid.to3d(grid.sp2grid(phi));

