function [phi, varargout] = groundstate_besp(task,dt,eps,phi0)
% groundstate_besp - Calculate the stationary state of GPE using Backwards
% Euler Spectral method
arguments
    task GPEtask
    dt double      %  evolution time step
    eps double     %  desired accuracy (applied to residual of the WF)
    phi0 = 'rand'  %  initial approximation of the wave function,
                   %  can be n-dim array of values or
                   %  'tf' - Thomas-Fermi initial approximation,
                   %  'rand' or empty - random
end
% Output argumants
%  phi   :  calculated stationary state
%  mu    :  (optional) history of chemical potential values from norm decrease
%  mu2   :  (optional) history of chemical potential from integral evaluation
%  E     :  (optional) history of energy values

grid = task.grid;
V = task.getVtotal(0);

phi = task.process_init_state(phi0);
task.current_state = phi;
 
MU = zeros(task.itp_max_iter,1,'like',V);
MU2 = zeros(task.itp_max_iter,1,'like',V);
EE = zeros(task.itp_max_iter,1,'like',V);
i = 0;
tmp = real(phi.*conj(phi)).*task.g+V;
while true
    i=i+1;
	bmax = max(tmp(:));
	bmin = min(tmp(:));
	alpha = (bmax+bmin)/2;
	phihat = grid.fft(phi);
	ghat = grid.fft((alpha-tmp).*phi);
	phi = grid.ifft((phihat+dt*ghat)./(1+dt*(alpha+grid.kk)));

    delta=1;
    while delta>eps
        ghat = grid.fft((alpha-tmp).*phi);
        phi11 = grid.ifft((phihat+dt*ghat)./(1+dt*(alpha+grid.kk)));
        delta = max(abs(abs(phi(:))-abs(phi11(:))));
        phi = phi11;
    end

    tmp = real(phi.*conj(phi));
    if(task.Ntotal > 0)
        mu = sqrt(task.Ntotal/grid.integrate(tmp));
        MU(i) = (mu-1)/dt;
    else
        mu = exp(task.mu_init*dt);
        MU(i) = grid.integrate(tmp*mu^2);
    end
    phi=phi*mu;
    tmp = tmp*mu^2.*task.g+V;
  
    if(i>task.itp_min_iter && mod(i,10) == 0)
        if(nargout >= 3)
            MU2(i) = real(grid.inner(phi,task.applyham(phi)));
        end        
        if(nargout >= 4)
            EE(i) = task.get_energy(phi);
        end  
        delta = max(abs(abs(phi(:))-abs(task.current_state(:))))/(dt*10);
        task.current_state = phi;
        if(delta < eps)
            break;
        end
    end
    
    if(i>=task.itp_max_iter)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:i);
    if(task.Ntotal > 0)
        task.current_mu = MU(end);
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    if(task.Ntotal > 0)
        MU2 = MU2(1:i)/task.Ntotal;
    else
        MU2 = MU2(1:i)./MU;
        task.current_mu = MU2(end);
        task.current_n = MU(end);
    end
    varargout{2} = MU2;
end
if(nargout >= 4)
    if(task.Ntotal > 0)
        EE = EE(1:i)./task.Ntotal;
    else
        EE = EE(1:i)./MU;
    end
    varargout{3} = EE;
end

task.init_state = phi;
end
