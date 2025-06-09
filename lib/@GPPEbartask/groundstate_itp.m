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
if task.bar_mass > 0
    task.bar_dens=task.particle_cloud_density()/task.grid.weight*task.bar_mass/task.bar_np;
    task.set_pot(task.bar_dens,true);
    tmp = task.Vtrap;
    task.Vtrap = task.Vtrap+task.Fi(task.grid_inds,task.grid_inds,task.grid_inds);
    task.set_pot(abs(phi0).^2,true);
end
[phi, varargout{1:nargout-1}] = task.groundstate_itp@GPPEtask(dt,eps,phi0);
if task.bar_mass > 0
    task.Vtrap = tmp;
end
end
