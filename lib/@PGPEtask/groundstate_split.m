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
g = task.g*task.Ntotal;
% omega = task.omega;
% n_cn=task.n_crank;
if(nargin <= 3)
    phi = grid.normalize(rand(size(grid.etot),'like',V) + 1i*rand(size(grid.etot),'like',V));
else
    phi = grid.normalize(phi0);
end
ekk = exp(-grid.etot*dt*0.5);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
delta = 1;
mu_old = 0;
i = 1;

while delta > eps
    
    phi = ekk.*phi;
    phir = grid.sp2grid(phi);
    phir = exp(-(V+g*abs(phir.^2))*dt).*phir;
    phi = grid.grid2sp(phir);
    phi = ekk.*phi;
    
    mu = sqrt(1.0/sum(abs(grid.to1d(phi.^2))));
    MU(i) = mu;
    
    phi=phi*mu;
    phir = grid.sp2grid(phi);

    if(nargout >= 3)
        MU2(i) = real(sum(grid.to1d(grid.etot.*abs(phi).^2 + conj(phi).*(grid.grid2sp((V+g*abs(phir.^2)).*phir)))));
    end
    if(i>50)
        delta = abs(log(mu_old/mu))/dt^2/10;
        mu_old = MU(i-10);
    end
    i=i+1;
    if(i>=10000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:nnz(MU));
    MU = 1/dt * log(MU);
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = MU2(1:nnz(MU2));
    varargout{2} = MU2;
end
phi = phi * sqrt(task.Ntotal);
task.init_state = phi;
phi=grid.to3d(grid.sp2grid(phi));

