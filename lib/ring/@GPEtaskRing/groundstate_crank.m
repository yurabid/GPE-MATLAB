function [phi, varargout] = groundstate_crank(task,dt,eps,phi0)
% groundstate_crank - Calculate the stationary state of GPE with Crank-Nicholson Imaginary Time Propagation method.
%
%  Usage :
%    phi = task.groundstate_crank(dt,eps)
%    phi = task.groundstate_crank(dt,eps,phi0)
%    [phi, mu] = task.groundstate_crank(dt,eps)
%    [phi, mu, mu2] = task.groundstate_crank(dt,eps)
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
omega = task.omega;
n_cn=task.n_crank;
if(nargin <= 3)
    phi = grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V));
else
    phi = grid.normalize(phi0);
end

MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
delta = 1;
mu_old = 0;
i = 1;
tmp2 = real(phi.*conj(phi));
while delta > eps && i<5000

    lphi = phi;
    for ii = 1:n_cn
        lphi = phi - dt*(task.applyh0(lphi,0) + g*tmp2.*lphi);
        lphi = 0.5*(phi+lphi);
    end
    phi = phi - dt*(task.applyh0(lphi,0) + g*tmp2.*lphi);

    
	tmp2 = real(phi.*conj(phi));
    mu = sqrt(1.0/grid.integrate(tmp2));
    phi=phi*mu;
	tmp2 = tmp2*mu^2;
    MU(i) = mu;
    if(nargout >= 3)
        MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + (V+g*tmp2).*tmp2));
        if(omega ~= 0)
            MU2(i) = MU2(i) - omega*real(grid.integrate(conj(phi).*grid.lz(phi)));
        end
    end
    if(i>50)
        delta = abs(log(mu_old/mu))/dt^2/10;
        mu_old = MU(i-10);
    end
    i=i+1;

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
task.current_mu = MU(end);
task.current_n = task.Ntotal;
end
