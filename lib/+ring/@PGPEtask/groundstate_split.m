function [phi, varargout] = groundstate_split(task,dt,eps,phi0)
% groundstate_split - Calculate the stationary state of GPE with split step Imaginary Time Propagation method.
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
% omega = task.omega;
% n_cn=task.n_crank;
if(nargin <= 3)
    phi = grid.normalize(rand(size(grid.etot),'like',V) + 1i*rand(size(grid.etot),'like',V));
else
    phi = grid.normalize(phi0);
end
ekk = exp(-grid.etot*dt*0.5);
ekkphi = exp(-grid.rdphimat*dt*0.5);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);
delta = 1;
mu_old = 0;
i = 0;
iswitch = 50;

while true
    
    i=i+1;
    phi = ekk.*phi;
    phir = grid.ifftr(phi);
    phir = ekkphi.*phir;
    phirzp = grid.ifftzphi(phir);
    phirzp = exp(-(V+g*abs(phirzp).^2)*dt).*phirzp;
    phir = grid.fftzphi(phirzp);
    phir = ekkphi.*phir;
    phi = grid.fftr(phir);
    phi = ekk.*phi;
    
%     mu = sqrt(1.0/sum(abs(grid.to1d(phi.^2))));
%     MU(i) = log(mu)/dt;

    if(task.mu_init > 0)
        mu = exp(task.mu_init*dt);
        MU(i) = sum(abs(grid.to1d(phi)).^2)*mu^2;
    else
        mu = sqrt(task.Ntotal/sum(abs(grid.to1d(phi)).^2));
        MU(i) = log(mu)/dt;
    end    
    
    phi=phi*mu;
    ncur = sum(abs(grid.to1d(phi)).^2);
%     phir = grid.sp2grid(phi);
    if(nargout >= 3)
        phir = grid.sp2grid(phi);
        h0 = real(sum(grid.to1d(conj(phi).*(grid.applyh0(phi) + grid.grid2sp(V.*phir)))));
%         h1 = real(sum(grid.to1d(conj(phi).*(grid.grid2sp(g*abs(phir).^2.*phir)))));
        h1 = grid.integrate_grid(g*abs(phir).^4);
        MU2(i) = (h0 + h1)/ncur;
        EE(i) = (h0 + h1*0.5)/ncur;
    end
%     if(i>50)
%         delta = abs(log(mu_old/mu))/dt^2/10;
%         mu_old = MU(i-10);
%     end
%     imagesc(real(phir(:,:,end/2))); drawnow;
    
    if((i-iswitch)>50 && mod(i,10) == 0)
%         delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt/MU(i);
        delta = abs((MU(i)-MU(i-10))^2/(MU(i)-2*MU(i-10)+MU(i-20)))/MU(i);
        if(delta < eps)
            if (dt<eps*10 || dt<1e-4)
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
%     MU = 1/dt * log(MU);
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
%phi = phi * sqrt(task.Ntotal);
task.init_state = phi;
phi=grid.to3d(grid.sp2grid(phi));

