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
%    phi0  :  initial approximation of the wave function,
%             'tf' - Thomas-Fermi initial approximation,
%             'rand' or empty - random
%  Output
%    phi      :  calculated stationary state
%    mu       :  array of chemical potential values from norm decrease
%    mu2      :  array of chemical potential from integral evaluation

grid = task.grid;
V = task.getVtotal(0);
g = task.g;
omega = task.omega;
n_cn=task.n_crank;
if(task.Ntotal > 0)
    nnn = task.Ntotal;
else
    nnn = 1;
end
if(nargin <= 3)
    phi0 = 'rand';
end
if(isa(phi0,'char'))
    if(task.Ntotal > 0)
        if(strcmp(phi0,'tf'))
            [phi,~] = task.groundstate_tf(eps); % Thomas-Fermi initial guess
        else
            phi = sqrt(nnn)*grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V)); % random initial guess
        end
    else
        phi = real(sqrt(complex(task.mu_init - V)/g)); % use only Thomas-Fermi approximation as initial guess if mu_init is set
    end
else
    phi = sqrt(nnn)*grid.normalize(phi0);
end

ekk = exp(-grid.kk*dt);
MU = zeros(5000,1,'like',V);
MU2 = zeros(5000,1,'like',V);
i = 0;

tmp2 = real(phi.*conj(phi))*g+V;
while true
    i=i+1;
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
    if(task.Ntotal > 0)
        mu = sqrt(task.Ntotal/grid.integrate(tmp));
        MU(i) = log(mu)/dt;
    else
        mu = exp(task.mu_init*dt);
        MU(i) = grid.integrate(tmp*mu^2);
    end
    phi=phi*mu;
    tmp = tmp*mu^2;
    tmp2 = tmp*g+V;
    
    if(nargout >= 3)
        MU2(i) = real(grid.integrate(conj(phi).*grid.lap(phi) + tmp.*tmp2));
        if(omega ~= 0)
            MU2(i) = MU2(i) - omega*real(grid.inner(phi,grid.lz(phi)));
        end
    end

    if(i>50 && mod(i,10) == 0)
        delta = (abs(MU(i)-MU(i-10))/10 + abs(MU(i)-MU(i-1)))/dt;
        if(delta < eps)
            break;
        end
    end
    
    if(i>=5000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:i);
    if(task.Ntotal > 0)
        MUEX = MU(i) - (MU(i)-MU(i-1))^2/(MU(i)-2*MU(i-1)+MU(i-2)); % exponential extrapolation
        MU = [MU; MUEX];
        task.current_mu = MUEX;
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    if(task.Ntotal > 0)
        MUEX = MU2(i) - (MU2(i)-MU2(i-1))^2/(MU2(i)-2*MU2(i-1)+MU2(i-2)); % exponential extrapolation
        MU2 = [MU2(1:i); MUEX]/task.Ntotal;
    else
        MU2 = MU2(1:i)./MU;
        task.current_mu = MU2(end);
        task.current_n = MU(end);
    end
    varargout{2} = MU2;
end

task.init_state = phi;
end
