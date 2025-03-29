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
if(nargin < 5)
    if(nargin < 4)
        if(nargin < 3)
            N = task.Ntotal;
        end
        V = task.getVtotal(0);
    end
    mul = min(V(:));
    mur = max(V(:));
end
NC = N*task.g(1);
mu = (mul+mur)/2;
while (mur-mul)/mu>eps
    vvv = (mu-V).*(V<mu);
    NN = task.grid.integrate(vvv);
    if(NN>NC)
        mur = mu;
    else
        mul = mu;
    end
    mu = (mul+mur)/2;
end
phi = sqrt(vvv/task.g(1));
