function [phi,mu] = groundstate_tf(task,eps,N,V,mul,mur)
% groundstate_tf - Calculate the Thomas-Fermi stationary state.
%
%  Usage :
%    [phi, mu] = task.groundstate_tf(eps)
%    [phi, mu] = task.groundstate_tf(eps,N,V,mul,mur)
%  Input
%    eps     :  desired accuracy (applied to chemical potential)
%    V       :  trap potential 
%    N       :  particle number
%    mur,mul :  upper and lower bounds on chemical potential
%  Output
%    phi      :  calculated stationary state
%    mu       :  chemical potential
grid = task.grid;

if(nargin < 5)
    if(nargin < 4)
        if(nargin < 3)
            N = task.Ntotal;
        end
        V = 0.5*(grid.omx^2*grid.mesh.x.^2+grid.omy^2*grid.mesh.y.^2+grid.omz^2*grid.mesh.z.^2) + task.getVtotal(0);
    end
    mul = min(V(:));
%     ndim = ndims(task.grid.mesh.x);
%     if(ndim==3)
    mur = min([V(1,end/2,end/2), V(end,end/2,end/2),V(end/2,1,end/2), V(end/2,end,end/2),V(end/2,end/2,1), V(end/2,end/2,end)]);
%     else
%         mur = min([V(1,end/2), V(end,end/2),V(end/2,1), V(end/2,end)]);
%     end
end

if(task.mu_init > 0)
    phi = real(sqrt(complex(task.mu_init - V)/task.g));
	mu = task.mu_init;
else
    grid2 = grid3d(grid.x(end),128,grid.y(end),128,grid.z(end),128);
    V2 = 0.5*(grid.omx^2*grid2.mesh.x.^2+grid.omy^2*grid2.mesh.y.^2+grid.omz^2*grid2.mesh.z.^2);
    NC = N*task.g/grid2.weight;
    mu = (mul+mur)/2;
    while (mur-mul)/mu>eps
        vvv = (mu-V2).*(V2<mu);
        NN = sum(vvv(:));
        if(NN>NC)
            mur = mu;
        else
            mul = mu;
        end
        mu = (mul+mur)/2;
    end
    phi = real(sqrt(complex(mu - V)/task.g));
%     phi = sqrt(vvv/task.g);
end
% phi=grid.grid2sp(phi);
