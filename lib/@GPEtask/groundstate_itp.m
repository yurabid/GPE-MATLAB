function [phi, varargout] = groundstate_itp(task,dt,eps)
% groundstate_itp - Calculate the stationary state of GPE with Imaginary Time Propagation method.
%
%  Usage :
%    phi = task.groundstate_itp(dt,eps)
%    [phi, mu] = task.groundstate_itp(dt,eps)
%    [phi, mu, mu2] = task.groundstate_itp(dt,eps)
%  Input
%    dt    :  evolution time step
%    eps   :  desired accuracy (applied to chemical potential)
%  Output
%    phi      :  calculated stationary state
%    mu       :  array of chemical potential values during evolution
%    mu2      :  array of chemical potential from integral evaluation

grid = task.grid;
V = task.getVtotal(0);
g = task.g*task.Ntotal;
omega = task.omega;
n_cn=10;

phi = grid.normalize(rand(size(grid.mesh.x),'like',grid.mesh.x));
ekk = exp(-grid.kk*dt);
MU = zeros(1000,1,'like',grid.mesh.x);
MU2 = zeros(1000,1,'like',grid.mesh.x);
delta = 1;
mu_old = 0;
i = 1;

tmp2 = abs(phi.*conj(phi));
while delta > eps
    phi = exp(-(V + g*tmp2)*dt*0.5).*phi;
    phi = grid.ifft(ekk.*grid.fft(phi));
        if(omega ~= 0)
            lphi = phi;
            for ii = 1:n_cn
                lphi = phi + dt*1i*omega*(grid.mesh.x.*grid.derivy(lphi) - ...
                    grid.mesh.y.*grid.derivx(lphi));
                lphi = 0.5*(phi+lphi);
            end
            phi = phi + dt*1i*omega*(grid.mesh.x.*grid.derivy(lphi) - ...
                grid.mesh.y.*grid.derivx(lphi));
        end
    phi = exp(-(V + g*phi.*conj(phi))*dt*0.5).*phi;
	tmp2 = abs(phi.*conj(phi));
    mu = sqrt(1.0/grid.integrate(tmp2));
    phi=phi*mu;
	tmp2 = tmp2*mu^2;
    MU(i) = mu;
    if(nargout >= 3)
        MU2(i) = real(grid.integrate(abs(conj(phi).*grid.ifft(grid.kk.*grid.fft(phi))) + (V+g*tmp2).*tmp2));
    end
    delta = abs(mu_old-mu)/dt;
    mu_old = mu;
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
end
