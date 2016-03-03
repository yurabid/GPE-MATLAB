function [phi, varargout] = groundstate_multi_itp(task,dt,eps,phi)
% groundstate_itp - Calculate the stationary state of GPE with Imaginary Time Propagation method.
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
V = task.getVtotal(0);
ncomp = size(task.g,1);
% nd = numel(size(grid.mesh.x));
g = task.g*task.Ntotal;
omega = task.omega;
n_cn=task.n_crank;
if(nargin <= 3)
    phi = cellfun(@(x) grid.normalize(rand(size(grid.mesh.x),'like',grid.mesh.x))./sqrt(ncomp), cell(1,ncomp), 'UniformOutput', false);
end
ekk = exp(-grid.kk*dt);
MU = zeros(1000,ncomp,'like',grid.mesh.x);
MU2 = zeros(1000,1,'like',grid.mesh.x);
delta = 1;
mu_old = 0;
i = 1;

tmp2 = num2cell(zeros(1,ncomp,'like',grid.mesh.x));
    for j =1:ncomp
        tmp2{j} = zeros(size(phi{j}),'like',grid.mesh.x);
        for k = 1:ncomp
            if(j~=k)
                tmp2{j} = tmp2{j} + g(j,k)*abs(phi{k}.*conj(phi{k}));
            end
        end
    end
    
while delta > eps && i<2000
    for j = 1:ncomp
        phi{j} = exp(-(V + tmp2{j}+g(j,j)*phi{j}.*conj(phi{j}))*dt*0.5).*phi{j};
        phi{j} = grid.ifft(ekk.*grid.fft(phi{j}));
            if(omega ~= 0)
                lphi = phi{j};
                for ii = 1:n_cn
                    lphi = phi + dt*1i*omega*(grid.mesh.x.*grid.derivy(lphi) - ...
                        grid.mesh.y.*grid.derivx(lphi));
                    lphi = 0.5*(phi+lphi);
                end
                phi{j} = phi{j} + dt*1i*omega*(grid.mesh.x.*grid.derivy(lphi) - ...
                    grid.mesh.y.*grid.derivx(lphi));
            end
        phi{j} = exp(-(V  + tmp2{j}+g(j,j)*phi{j}.*conj(phi{j}))*dt*0.5).*phi{j};
    end

    if(ncomp > 1)
        for j =2:ncomp
            curnorm = grid.norm(phi{j});
            for k = 1:j-1
                phi{j} = phi{j} - grid.inner(phi{k},phi{j})/grid.inner(phi{k},phi{k}) * phi{k};
            end
            phi{j} = phi{j}*sqrt(curnorm/grid.norm(phi{j}));
        end
    end
    
    dens = zeros(size(phi{j}),'like',grid.mesh.x);
    for j =1:ncomp
        tmp2{j} = zeros(size(phi{j}),'like',grid.mesh.x);
        for k = 1:ncomp
            if(j~=k)
                tmp2{j} = tmp2{j} + g(j,k)*abs(phi{k}.*conj(phi{k}));
            end
        end
        dens = dens + abs(phi{k}.*conj(phi{k}));
    end
    
%     mu = sqrt(1.0/grid.integrate(dens));
    mu = zeros(1,ncomp,'like',grid.mesh.x);
    for j =1:ncomp
        mu(j) = sqrt(1.0/grid.integrate(abs(phi{j}.*conj(phi{j})))/ncomp);
        phi{j}=phi{j}*mu(j);
        tmp2{j} = tmp2{j}*mu(j)^2;
    end
    MU(i,:) = mu;
    

    
    if(nargout >= 3)
%         MU2(i) = real(grid.integrate(abs(conj(phi).*grid.ifft(grid.kk.*grid.fft(phi))) + (V+g*tmp2).*tmp2));
    end
    if(i>50)
        delta = sum(abs(log(mu_old(:)./mu(:)))/dt^2*10);
        mu_old = MU(i-10,:);
    end
    i=i+1;
end

if(nargout >= 2)
    MU = MU(1:nnz(MU(:,1)),:);
    MU = 1/dt * log(MU);
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = MU2(1:nnz(MU2));
    varargout{2} = MU2;
end
for j =1:ncomp
    phi{j} = phi{j} * sqrt(task.Ntotal);
end
task.init_state = phi;
end
