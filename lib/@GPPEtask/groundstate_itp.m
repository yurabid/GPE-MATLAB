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
h=grid.x(2)-grid.x(1);
sz = grid.nx+1;
if(nargin <= 3)
    phi0 = 'rand';
end
if(isa(phi0,'char'))
    phi = sqrt(task.Ntotal)*grid.normalize(rand(size(grid.mesh.x),'like',grid.mesh.x) + 1i*rand(size(grid.mesh.x),'like',grid.mesh.x)); % random initial guess
else
    phi = sqrt(task.Ntotal)*grid.normalize(phi0);
end
if(numel(task.Fi)==0)
    task.set_pot(phi);
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
            [task.Fi,~] = task.V_cycle(task.Fi,task.nlinop,h,sz);
        end
    end
    task.current_iter = i;

    if(i>50 && mod(i,10) == 5)
        if(nargout >= 3)
            MU2(i) = real(grid.inner(phi,task.applyham(phi)));
        end
        if(nargout >= 4)
            EE(i) = task.get_energy(phi)/task.Ntotal;
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
    MU = MU(1:i);
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = MU2(1:i)/task.Ntotal;
    varargout{2} = MU2;
end
if(nargout >= 4)
    EE = EE(1:i);
    varargout{3} = EE;
end
end
