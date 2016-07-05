function [phi, varargout] = groundstate_itp(task,dt,eps,phi0)
% groundstate_itp - Calculate the stationary state of GPE with split step Imaginary Time Propagation method.
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
g = task.g*task.Ntotal;
omega = task.omega;
n_cn=task.n_crank;
if(nargin <= 3)
    phi = grid.normalize(rand(size(grid.mesh.x),'like',grid.mesh.x) + 1i*rand(size(grid.mesh.x),'like',grid.mesh.x));
else
    phi = grid.normalize(phi0);
end
ekk = exp(-grid.kk*dt);
MU = zeros(1000,1,'like',grid.mesh.x);
MU2 = zeros(1000,1,'like',grid.mesh.x);
delta = 1;
mu_old = 0;
i = 1;

tmp2 = real(phi.*conj(phi))*g+V;
while delta > eps
    phi = exp(-tmp2*dt*0.5).*phi;
    phi = grid.ifft(ekk.*grid.fft(phi));
        if(omega ~= 0)
            lphi = phi;
            for ii = 1:n_cn
                lphi = phi + dt*omega*grid.lz(lphi);
                lphi = 0.5*(phi+lphi);
            end
            phi = phi + dt*omega*grid.lz(lphi);
        end
    phi = exp(-tmp2*dt*0.5).*phi;
    
	tmp = real(phi.*conj(phi));
    mu = sqrt(1.0/grid.integrate(tmp));
    phi=phi*mu;
	tmp = tmp*mu^2;
	tmp2 = tmp*g+V;
    MU(i) = mu;
    if(nargout >= 3)
        MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + tmp.*tmp2));
        if(omega ~= 0)
            MU2(i) = MU2(i) - omega*real(grid.integrate(conj(phi).*grid.lz(phi)));
        end
    end
    if(i>50)
        delta = abs(log(mu_old/mu))/dt^2/10;
        mu_old = MU(i-10);
    end
    i=i+1;
    if(i>=5000) 
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
end
