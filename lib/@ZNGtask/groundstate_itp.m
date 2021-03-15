function [phi, varargout] = groundstate_itp(task,dt,eps,phi0)
% groundstate_itp - Calculate the Hartree-Fock stationary state with split step Imaginary Time Propagation method.
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

TT = task.T;

grid = task.grid;
grid_nt = task.grid_nt;
nx = grid.nx/2;
ny = grid.ny/2;
nz = grid.nz/2;
V = task.getVtotal(0);
g = task.g;
omega = task.omega;
n_cn=task.n_crank;

if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = 1;
end
if(nargin <= 3)
    phi0 = 'rand';
end
if(isa(phi0,'char'))
    if(task.Ntotal > 0)
        if(strcmp(phi0,'tf'))
            phi = task.groundstate_tf(eps,1); % Thomas-Fermi initial guess
        else
            phi = sqrt(NN0)*grid.normalize(rand(size(grid.mesh.x),'like',V) + 1i*rand(size(grid.mesh.x),'like',V)); % random initial guess
        end
    else
        phi = real(sqrt(complex(task.mu_init - V)/g)); % use only Thomas-Fermi approximation as initial guess if mu_init is set
    end   
else
    phi = sqrt(NN0)*grid.normalize(phi0);
end
nt = zeros(size(grid_nt.mesh.x),'like',V);
tmpext = task.vtrap_nt;%zeros(size(grid_nt.mesh.x),'like',V);
NNN = NN0;
NNt = 0;
ekk = exp(-grid.kk*dt);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
delta = 1;
mu_old = 0;
i = 1;

tmp2 = real(phi.*conj(phi))*g+V;
while (delta > eps || mod(i,10)~=5)
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
        NNt = grid_nt.integrate(nt);
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
    MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + tmp.*(tmp2+2*g*nt(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz))));
    if(omega ~= 0)
        MU2(i) = MU2(i) - omega*real(grid.integrate(conj(phi).*grid.lz(phi)));
    end
    MU2(i) = MU2(i)/ncur;
    if(mod(i,10)==0)%i>150)
		if(TT>0)
		    tmpext(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz) = V + 2*g*tmp;
            mmu = min(MU2(i),min(min(min(tmpext + 2*g*nt)))-1e-10); % compensate for possibly inaccurate chem. pot. calculation
            
		    ntt=(TT/(2*pi))^1.5*polylog(1.5,exp((mmu - tmpext - 2*g*nt)/TT)); % averaging increases stability for high temperatures
		    if(NNN<=1 && task.Ntotal > 0)
		        NNt = grid_nt.integrate(ntt);
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
	tmp2 = tmp2 + 2*g*nt(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz);
    i=i+1;
    if(i>=10000) 
        warning('Convergence not reached');
        break;
    end
end
i=i-1;
if(nargout >= 2)
    MU = real(MU(1:nnz(MU)));
    if(task.Ntotal > 0)
        MU = 1/dt * log(MU);
        task.current_mu = MU(end);
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = real(MU2(1:nnz(MU2)));
%     MUEX = MU2(i) - (MU2(i)-MU2(i-1))^2/(MU2(i)-2*MU2(i-1)+MU2(i-2)); % exponential extrapolation
    varargout{2} = MU2;
end

task.init_state = phi;
task.init_state_nt = nt;

end
