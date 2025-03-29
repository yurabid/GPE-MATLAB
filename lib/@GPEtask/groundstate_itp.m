function [phi, varargout] = groundstate_itp(task,dt,eps,phi0)
% groundstate_itp - Calculate the stationary state of GPE with split-step 
% Imaginary Time Propagation method.
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
if(task.Ntotal > 0)
    nnn = task.Ntotal;
else
    nnn = 1;
end

phi = task.process_init_state(phi0);
task.current_state = phi;

task.set_kinop(dt);
MU = zeros(task.itp_max_iter,1,'like',V);
MU2 = zeros(task.itp_max_iter,1,'like',V);
EE = zeros(task.itp_max_iter,1,'like',V);
i = 0;
dtmin = sqrt(eps);
task.nlinop = exp(-dt*0.5*(abs(phi).^2.*task.g+V));
while true
    i=i+1;

    phi = task.nlinop.*phi;
    phi = task.ssft_kin_step(phi,dt);
    phi = task.nlinop.*phi;

    task.nlinop = abs(phi).^2;
    if(task.Ntotal > 0)
        mu = sqrt(task.Ntotal/grid.integrate(task.nlinop));
        MU(i) = log(mu)/dt;
    else
        mu = exp(task.mu_init*dt);
        MU(i) = grid.integrate(task.nlinop)*mu^2;
    end
    phi=phi*mu;

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
            if(nargout >= 3)
                dtmin = dt*sqrt(eps/abs(MU(i)-MU2(i)/nnn));
            end
            if (~task.itp_adjust_stepsize || dt<dtmin)
                break;
            else
                dt = dt/2;
                task.set_kinop(dt);
            end
        end
    end
    
    task.nlinop = exp(-dt*0.5*(task.nlinop.*task.g*mu^2+V));

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
