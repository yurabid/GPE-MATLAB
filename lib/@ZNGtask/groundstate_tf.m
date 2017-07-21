function [phi, varargout] = groundstate_tf(task,eps,dummy)
% groundstate_itp - Calculate the Thomas-Fermi stationary state.
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
mu_min = min(V(:));
% mu_max = min([V(1,end/2,end/2), V(end,end/2,end/2),V(end/2,1,end/2), V(end/2,end,end/2),V(end/2,end/2,1), V(end/2,end/2,end)]);
mu_max = max(V(:));
g = task.g;
TT = task.T;
NN0 = task.Ntotal;
[phi,mu] = task.groundstate_tf@GPEtask(eps,NN0,V,mu_min,mu_max);
if(nargin == 3)
    return
end
nt = phi*0;
MU = zeros(1000,1,'like',grid.mesh.x);
delta = 1;
mu_old = 0;
i = 1;
while delta > eps
    mmu = min(mu,min(min(min(V+2*g*(phi.^2+nt))))-1e-10); % make sure that polylog will not brake
    ntt=(TT/(2*pi))^1.5*polylog(1.5,exp((mmu-V-2*g*(phi.^2+nt))/TT));
    nt = (nt+ntt)/2;
    NNt = grid.integrate(nt);
    NNN = NN0 - NNt;
    if(NNN<1)
        NNN=1;
        nt = nt*NN0/NNt; % we need to get the correct total number of particles even above Tc
    end
    [phi,mu] = task.groundstate_tf@GPEtask(eps,NNN,V+2*g*nt,mu_min,mu_max);
    MU(i)=mu;
    if(i>20)
        delta = abs(MU(i)-mu_old)/10;
        mu_old = MU(i-10);
    end
    i=i+1;
    if(i>=5000) 
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = real(MU(1:nnz(MU)));
    varargout{1} = MU;
end

task.init_state = phi;
task.init_state_nt = nt;
task.current_mu = MU(end);
task.current_n = task.Ntotal;
