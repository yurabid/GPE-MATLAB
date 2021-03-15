function [phi, varargout] = groundstate(task,eps,phi0)
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
rscale = grid.weight^(1/grid.ndims);
escale = 1/rscale^2;
% tscale = rscale^2;
phiscale = rscale^(-grid.ndims/2);
kk = grid.kk/escale;
V = task.getVtotal(0)/escale;
g = task.g*phiscale^2/escale;
dt = 0.5;
omega = task.omega/escale;
% n_cn=task.n_crank;
if(task.Ntotal > 0)
    nnn = task.Ntotal;
else
    nnn = 1;
end
if(nargin <= 2)
    phi0 = 'rand';
end
if(isa(phi0,'char'))
    if(task.Ntotal > 0)
        if(strcmp(phi0,'tf'))
            [phi,~] = task.groundstate_tf(eps); % Thomas-Fermi initial guess
        else
            phi = sqrt(nnn)*grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V)); % random initial guess
        end
    else
        phi = real(sqrt(complex(task.mu_init - V)./g)); % use only Thomas-Fermi approximation as initial guess if mu_init is set
    end
else
    phi = sqrt(nnn)*grid.normalize(phi0);
end

phi = phi/phiscale;

ekk = exp(-kk*0.5*dt);
if(omega ~= 0)
    ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*task.omega)/escale/4*dt);
    eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*task.omega)/escale/4*dt);
end
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);
% dts=[];
% deltas = [];
i = 0;
iswitch = 50;

tmp2 = real(phi.*conj(phi)).*g+V;
while true
    i=i+1;
%     phi = exp(-tmp2*dt*0.5).*phi;
    if(omega ~= 0)
        phi = grid.ifftx(ekx.*grid.fftx(phi));
        phi = grid.iffty(eky.*grid.ffty(phi));
    else
        phi = grid.ifft(ekk.*grid.fft(phi));
    end
        phi = exp(-tmp2*dt).*phi;

    if(omega ~= 0)
        phi = grid.iffty(eky.*grid.ffty(phi));
        phi = grid.ifftx(ekx.*grid.fftx(phi));
    else
        phi = grid.ifft(ekk.*grid.fft(phi));
    end
%     phi = exp(-tmp2*dt*0.5).*phi;
    
    tmp = real(phi.*conj(phi));
    if(task.Ntotal > 0)
        mu = sqrt(task.Ntotal/grid.integrate(tmp))/phiscale;
        MU(i) = log(mu)/dt*escale;
    else
        mu = exp(task.mu_init/escale*dt);
        MU(i) = grid.integrate(tmp)*mu^2*phiscale^2;
    end
    phi=phi*mu;
    tmp = tmp*mu^2;
% imagesc(abs(phi));drawnow;
    tmp2 = tmp.*g+V;

    if(nargout >= 3)
        ncur = 1;%grid.integrate(tmp)*phiscale^2;
        hphi = task.applyham(phi*phiscale);
        MU2(i) = real(grid.inner(phi*phiscale,hphi))/ncur;
    else
        MU2(i) = MU(i);
    end
    if(nargout >= 4)
        EE(i) = task.get_energy(phi*phiscale)/ncur;
%         EE(i) = grid.integrate(abs(MU2(i)*phi*phiscale-hphi).^2)/task.Ntotal; 
    else
        EE(i) = MU2(i);
    end
% dts(i) = dt;
% if(i>20)
% deltas(i) = abs((EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)))/EE(i);
% semilogy(deltas);drawnow;
% end
    if((i-iswitch)>100 && mod(i,10) == 0)
        delta = abs((EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)))/EE(i);
%         if((EE(i)-2*EE(i-10)+EE(i-20))<0)
%             phi = phi + 0.1*grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V));
%         end
%         delta = abs((MU2(i)-MU2(i-10))^2/(MU2(i)-2*MU2(i-10)+MU2(i-20)))/MU2(i);
%         delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt/MU(i);
%         delta = EE(i);
        if(delta < eps)
            if (dt<eps*10 || dt<1e-4)
                break;
            else
                dt = dt/1.5;
                ekk = exp(-0.5*kk*dt);
                if(omega ~= 0)
                    ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*task.omega)/escale/4*dt);
                    eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*task.omega)/escale/4*dt);
                end
                iswitch = i;
            end
        end
    end
    
%     semilogy(EE);hold on;plot((MU2(:)-MU(:))./dts(:).^2);hold off;drawnow;
    
    if(i>=10000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:i);
    if(task.Ntotal > 0)
%         MUEX = MU(i) - (MU(i)-MU(i-10))^2/(MU(i)-2*MU(i-10)+MU(i-20)); % exponential extrapolation
%         MU = [MU; MUEX];
        task.current_mu = MU(end);
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    if(task.Ntotal > 0)
        MU2 = MU2(1:i);
%         MUEX = MU2(i) - (MU2(i)-MU2(i-10))^2/(MU2(i)-2*MU2(i-10)+MU2(i-20)); % exponential extrapolation
%         MU2 = [MU2(1:i); MUEX];
    else
        MU2 = MU2(1:i)./MU;
        task.current_mu = MU2(end);
        task.current_n = MU(end);
    end
    varargout{2} = MU2;
end

if(nargout >= 4)
    if(task.Ntotal > 0)
        EE = EE(1:i);
%         MUEX = EE(i) - (EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)); % exponential extrapolation
%         EE = [EE(1:i); MUEX];
    else
        EE = EE(1:i)./MU;
    end
    varargout{3} = EE;
end

phi = phi*phiscale;
task.init_state = phi;
task.current_state = phi;
end
