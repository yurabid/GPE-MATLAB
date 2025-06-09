function [phi, varargout] = groundstate_itp(task,dt,eps,phi0)
% groundstate_itp - Calculate the stationary state of GPE with split-step 
% Imaginary Time Propagation method.
arguments
    task GPEtask
    dt double      %  evolution time step
    eps double     %  desired accuracy (applied to residual of the WF)
    phi0 = 'rand'  %  initial approximation of the wave function,
                   %  can be n-dim array of values or
                   %  'rand' or empty - random
end
% Output argumants
%  phi   :  calculated stationary state
%  mu    :  (optional) history of chemical potential values from norm decrease
%  mu2   :  (optional) history of chemical potential from integral evaluation
%  E     :  (optional) history of energy values

grid = task.grid;
sz = grid.nx+1;
phi = task.process_init_state(phi0);
if(numel(task.Fi)==0)
    task.set_pot(abs(phi).^2);
end
task.current_state = phi;
if(task.itp_use_fft_Fi)
    kk_inv = 1./grid.kk/2;
    kk_inv(1,1,1) = 0;
end
  
task.set_kinop(dt);
MU = zeros(task.itp_max_iter,1,'like',grid.mesh.x);
MU2 = zeros(task.itp_max_iter,1,'like',grid.mesh.x);
EE = zeros(task.itp_max_iter,1,'like',grid.mesh.x);
i = 0;

task.nlinop = exp(-dt*0.5*(real(phi.*conj(phi)).*task.g+task.getVtotal(0)));
while true
    i=i+1;
    phi = task.nlinop.*phi; 
    phi = task.ssft_kin_step(phi,dt);
    phi = task.nlinop.*phi;

    task.nlinop = real(phi.*conj(phi));
    mu = sqrt(task.Ntotal/grid.integrate(task.nlinop));
    MU(i) = log(mu)/dt;
    phi=phi*mu;
    task.nlinop = task.nlinop*mu^2;
    
    if(mod(i,task.itp_fi_update_step) == 0)
        if(task.itp_use_fft_Fi)
            task.Fi = -grid.ifft(kk_inv.*grid.fft(task.nlinop));
        else
            task.set_pot_bc(task.nlinop);
            [task.Fi,~] = task.V_cycle(task.Fi,task.nlinop,grid.dx,sz);
        end
    end
    task.current_iter = i;

    if(i>task.itp_min_iter && mod(i,10) == 5)
        if(nargout >= 3)
            MU2(i) = real(grid.inner(phi,task.applyham(phi)));
        end
        if(nargout >= 4)
            EE(i) = task.get_energy(phi);
        end        

        delta = max(abs(abs(phi(:))-abs(task.current_state(:))))/(dt*10);
        task.current_state = phi;
        if(delta < eps)
            if (~task.itp_adjust_stepsize || dt<sqrt(eps))
                break;
            else
                dt = dt/2;
                task.set_kinop(dt);
            end
        end
    end
    
    task.nlinop = exp(-dt*0.5*(task.nlinop.*task.g+task.getVtotal(0)));    
    
    if(i>=task.itp_max_iter)
        warning('Convergence not reached');
        break;
    end
end

task.current_mu = MU(i);
task.current_n = task.Ntotal;
    
if(nargout >= 2)
    varargout{1} = MU(1:i);
end
if(nargout >= 3)
    varargout{2} = MU2(1:i)/task.Ntotal;
end
if(nargout >= 4)
    varargout{3} = EE(1:i)/task.Ntotal;
end
end
