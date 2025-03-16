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
            % Thomas-Fermi initial guess
            [phi,~] = task.groundstate_tf(eps); 
        else
            % random initial guess
            phi = sqrt(nnn)*grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V)); 
        end
    else
        % use only Thomas-Fermi approximation as initial guess if mu_init is set
        phi = real(sqrt(complex(task.mu_init - V)./task.g));
    end
else
    phi = sqrt(nnn)*grid.normalize(phi0);
end
task.current_state = phi;

task.set_kinop(dt);
MU = zeros(task.itp_max_iter,1,'like',V);
MU2 = zeros(task.itp_max_iter,1,'like',V);
EE = zeros(task.itp_max_iter,1,'like',V);
i = 0;
dtmin = 0;
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

    if(i>50 && mod(i,10) == 0)
        if(nargout >= 3)
            MU2(i) = real(grid.inner(phi,task.applyham(phi)));
        end
        if(nargout >= 4)
            EE(i) = task.get_energy(phi);
        end
        % delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt;
        delta = max(abs(abs(phi(:))-abs(task.current_state(:))))/(dt*10);
        task.current_state = phi;
        if(delta < eps)
            if(dtmin == 0)
                if(nargout >= 3)
                    dtmin = dt*sqrt(eps/abs(MU(i)-MU2(i)/nnn));
                else
                    dtmin = sqrt(eps);
                end
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
