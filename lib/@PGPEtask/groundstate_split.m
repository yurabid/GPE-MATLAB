function [phi, varargout] = groundstate_split(task,dt,eps,phi0)
% groundstate_split - Calculate the stationary state of GPE with split step Imaginary Time Propagation method.
%
%  Usage :
%    phi = task.groundstate_itp(dt,eps)
%    phi = task.groundstate_itp(dt,eps,phi0)
%    [phi, mu] = task.groundstate_itp(dt,eps)
%    [phi, mu, mu2] = task.groundstate_itp(dt,eps)
%  Input
%    dt    :  evolution time step
%    eps   :  desired accuracy (applied to chemical potential)
%    phi0  :  initial approximation of the wave function
%  Output
%    phi      :  calculated stationary state
%    mu       :  array of chemical potential values from norm decrease
%    mu2      :  array of chemical potential from integral evaluation

grid = task.grid;
V = task.getVtotal(0);
task.V0 = V;
g = task.g*task.Ntotal;
if(nargin <= 3)
    phi = grid.normalize(rand(size(grid.etot),'like',V) + 1i*rand(size(grid.etot),'like',V));
else
    phi = grid.normalize(phi0);
end
ekk = exp(-grid.etot*dt*0.5);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);

i = 0;
phir = grid.sp2grid(phi);
tmp = abs(phir.^2);

while true %delta > eps
    i=i+1;
   
    phi = ekk.*phi;
    phir = grid.sp2grid(phi);
    phir = exp(-(V+g*tmp)*dt).*phir;
    phi = grid.grid2sp(phir);
    phi = ekk.*phi;

    mu = sqrt(1.0/sum(abs(grid.to1d(phi.^2))));
    MU(i) = log(mu)/dt;
    
    phi=phi*mu;
    phir = grid.sp2grid(phi);
    tmp = abs(phir.^2);

    if(nargout >= 3)
        h1 = real(sum(grid.to1d(grid.etot.*abs(phi).^2 + conj(phi).*(grid.grid2sp(V.*phir)))));
        h2 = real(sum(conj(phi).*(grid.grid2sp(g*tmp.*phir))));
        MU2(i) = h1+h2;
        EE(i) = h1+h2/2;
    end
    if((i)>50 && mod(i,10) == 0)
%         delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt/MU(i);
        delta = abs((MU(i)-MU(i-10))^2/(MU(i)-2*MU(i-10)+MU(i-20)))/MU(i);
        if(delta < eps)
            if (dt<eps || dt<1e-5)
                break;
            else
                dt = dt/1.2;
                ekk = exp(-grid.etot*dt*0.5);
            end
        end
    end
    if(i>=10000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:nnz(MU));
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = MU2(1:nnz(MU2));
    varargout{2} = MU2;
end
if(nargout >= 4)
    EE = EE(1:nnz(EE));
    varargout{3} = EE;
end
phi = phi * sqrt(task.Ntotal);
task.init_state = phi;
phi=grid.to3d(grid.sp2grid(phi));

