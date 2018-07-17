function [phi, varargout] = groundstate_itp(task,dt,eps,phi)
% groundstate_itp - Calculate the stationary state of multicomponent GPE with split-step Imaginary Time Propagation method.
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
%    mu       :  array of chemical potential values during evolution
%    mu2      :  array of chemical potential from integral evaluation

grid = task.grid;
ncomp = size(task.g,1);
task.ncomp = ncomp;
VV = task.getVtotal(0);
V = VV{1};
coupl = task.coupling;
g = task.g;
if(nargin <= 3)
    phi = cellfun(@(x) grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V))./sqrt(ncomp), cell(1,ncomp), 'UniformOutput', false);
end
ekk = exp(-grid.kk*dt*0.5);
MU = zeros(1000,'like',V);
delta = 1;
i = 1;

tmp = cell(1,ncomp);
tmp2 = cell(1,ncomp);
% tmp3 = cell(1,ncomp);
cang = angle(coupl);
cosom = cosh(dt*abs(coupl));
sinomm = sinh(dt*abs(coupl)).*exp(-1i*cang);
sinomp = sinh(dt*abs(coupl)).*exp(1i*cang);

    
while delta > eps

    for j =1:ncomp
        tmp2{j} = VV{j};
        for k = 1:ncomp
            tmp2{j} = tmp2{j} + g(j,k)*abs(phi{k}).^2;%.*conj(phi{k}));
        end
    end     
    for j = 1:ncomp
        phi{j} = grid.ifft(ekk.*grid.fft(phi{j}));
    end
    %for j = 1:ncomp
        tmp2{1} = exp(-tmp2{1}*dt).*(cosom.*phi{1} - sinomm.*phi{2});
        tmp2{2} = exp(-tmp2{2}*dt).*(cosom.*phi{2} - sinomp.*phi{1});
    %end

    for j = 1:ncomp
        phi{j} = grid.ifft(ekk.*grid.fft(tmp2{j}));
    end
    
    ntot = 0;
    for j =1:ncomp
        tmp{j} = real(phi{j}.*conj(phi{j}));
        ntot = ntot + grid.integrate(tmp{j});
    end
    mu = sqrt(task.Ntotal/ntot);
    MU(i) = mu;
    
    for j =1:ncomp
        phi{j}=phi{j}.*mu;
    end
    
    if(i>50)
        delta = sum(abs(log(MU(i-10)./mu))/dt^2*10);
    end
    i=i+1;
    if(i>=10000) 
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:nnz(MU));
    MU = log(MU)/dt;
    varargout{1} = MU;
end
if(nargout >= 3)
    mu = 0;
    hphi = task.applyham(phi);
    for i =1:ncomp
        mu = mu + real(grid.inner(phi{i},hphi{i}));
    end
    varargout{2} = mu/task.Ntotal;    
end
task.init_state = phi;
task.current_state = phi;
end
