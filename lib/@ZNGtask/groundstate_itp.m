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
g = task.g;
omega = task.omega;
n_cn=task.n_crank;
TT = task.T;
if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = 1;
end
if(nargin <= 3)
    phi = sqrt(NN0)*grid.normalize(rand(size(grid.mesh.x),'like',grid.mesh.x) + 1i*rand(size(grid.mesh.x),'like',grid.mesh.x));
else
    phi = sqrt(NN0)*grid.normalize(phi0);
end
nt = phi*0;
NNN = NN0;
NNt = 0;
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
    if(TT>0 && mod(i,10)==1 && task.Ntotal > 0) % for better performance and stability we do some initial iterations without a thermal cloud
        NNt = grid.integrate(nt);
        NNN = NN0 - NNt;
        if(NNN<1)
            NNN=1;
            nt = nt*NN0/NNt; % we need to get the correct total number of particles even above Tc
        end
    end
    ncur = grid.integrate(tmp);
    if(task.Ntotal > 0)
        mu = sqrt(NNN/ncur);
        MU(i) = mu;
        ncur = ncur*mu^2;
    else
        mu = exp(task.mu_init*dt);
        ncur = ncur*mu^2;
        MU(i) = ncur;
    end
    phi=phi*mu;
	tmp = tmp*mu^2;
	tmp2 = tmp*g+V;
    MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + tmp.*(tmp2+2*g*nt)));
    if(omega ~= 0)
        MU2(i) = MU2(i) - omega*real(grid.integrate(conj(phi).*grid.lz(phi)));
    end
    MU2(i) = MU2(i)/ncur;
    if(mod(i,10)==0)%i>150)
		if(TT>0)
		    mmu = min(MU2(i),min(min(min(V+2*g*(tmp+nt))))-1e-10); % compensate for possibly inaccurate chem. pot. calculation
		    ntt=(TT/(2*pi))^1.5*polylog(1.5,exp((mmu-V-2*g*(tmp+nt))/TT)); % averaging increases stability for high temperatures
		    if(NNN<=1 && task.Ntotal > 0)
		        NNt = grid.integrate(ntt);
		        nt = (nt+ntt*NN0/NNt)*0.5;
		    else
		        nt=(nt+ntt)*0.5;
		    end
		end
        if(task.Ntotal > 0)
            delta = abs(log(mu_old/mu))/dt^2/9;
        else
            delta = abs(MU(i)-mu_old)/dt/9;
        end
        mu_old = MU(i-9);
    end
	tmp2 = tmp2 + 2*g*nt;
    i=i+1;
    if(i>=10000) 
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = real(MU(1:nnz(MU)));
    if(task.Ntotal > 0)
        MU = 1/dt * log(MU);
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = real(MU2(1:nnz(MU2)));
    varargout{2} = MU2;
end

task.init_state = phi;
task.init_state_nt = nt;
end
