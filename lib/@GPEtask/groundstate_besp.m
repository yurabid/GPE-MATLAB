function [phi, varargout] = groundstate_besp(task,dt,eps,phi0)
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
% n_cn=task.n_crank;
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
        phi = real(sqrt(complex(task.mu_init - V)./g)); % use only Thomas-Fermi approximation as initial guess if mu_init is set
    end
else
    phi = sqrt(nnn)*grid.normalize(phi0);
end

%ekk = exp(-grid.kk*dt);
%if(omega ~= 0)
%    ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*task.omega)/4*dt);
%    eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*task.omega)/4*dt);
%end    
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);
i = 0;
iswitch = 20;
eps2=eps*100;
tmp2 = real(phi.*conj(phi)).*g+V;
while true
    i=i+1;
	bmax = max(tmp2(:));
	bmin = min(tmp2(:));
	alpha = (bmax+bmin)/2;

	phihat = grid.fft(phi);
	ghat = grid.fft((alpha-tmp2).*phi);
	phi1 = grid.ifft((phihat+dt*ghat)./(1+dt*(alpha+grid.kk)));

    delta=1;
    while delta>eps2
        ghat = grid.fft((alpha-tmp2).*phi1);
        phi11 = grid.ifft((phihat+dt*ghat)./(1+dt*(alpha+grid.kk)));
        delta = max(abs(phi1(:)-phi11(:)));
        phi1 = phi11;
    end

    
    tmp = real(phi1.*conj(phi1));
    if(task.Ntotal > 0)
        mu = sqrt(task.Ntotal/grid.integrate(tmp));
        MU(i) = (mu-1)/dt;
    else
        mu = exp(task.mu_init*dt);
        MU(i) = grid.integrate(tmp*mu^2);
    end
    phi1=phi1*mu;
    tmp = tmp*mu^2;
    tmp2 = tmp.*g+V;

    if(nargout >= 3)
        MU2(i) = real(grid.inner(phi1,task.applyham(phi1)));
    else
        MU2(i) = MU(i);
    end        

    if(nargout >= 4)
        EE(i) = task.get_energy(phi1)/task.Ntotal;
    else
        EE(i) = MU2(i);
    end    
    if i==iswitch
        eps2=eps;
    end
    if(i>iswitch && mod(i,10) == 0)
%         delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt/MU(i);
%         delta = abs((EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)))/EE(i);
        delta = max(abs(phi(:)-phi1(:)));
        if(delta < eps)
            break;
        end
    end
    phi=phi1;
    if(i>=50000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:i);
    if(task.Ntotal > 0)
        MUEX = MU(i) - (MU(i)-MU(i-5))^2/(MU(i)-2*MU(i-5)+MU(i-10)); % exponential extrapolation
        MU = [MU; MUEX];
        task.current_mu = MUEX;
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    if(task.Ntotal > 0)
        MUEX = MU2(i) - (MU2(i)-MU2(i-5))^2/(MU2(i)-2*MU2(i-5)+MU2(i-10)); % exponential extrapolation
        MU2 = [MU2(1:i); MUEX]/task.Ntotal;
    else
        MU2 = MU2(1:i)./MU;
        task.current_mu = MU2(end);
        task.current_n = MU(end);
    end
    varargout{2} = MU2;
end
if(nargout >= 4)
    if(task.Ntotal > 0)
        MUEX = EE(i) - (EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)); % exponential extrapolation
        EE = [EE(1:i); MUEX];
    else
        EE = EE(1:i)./MU;
    end
    varargout{3} = EE;
end

task.init_state = phi;
task.current_state = phi;
end
