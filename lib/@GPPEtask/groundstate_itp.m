function [phi, varargout] = groundstate_itp(task,dt,eps,phi0,tc)
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
g = task.g;
h=grid.x(2)-grid.x(1);
angs=atan2(grid.mesh.y,grid.mesh.x);
sz = size(grid.mesh.x)+1;
omega = task.omega;
nnn = task.Ntotal;
if(nargin <= 3)
    phi0 = 'rand';
end
if(isa(phi0,'char'))
    phi = sqrt(nnn)*grid.normalize(rand(size(grid.mesh.x),'like',grid.mesh.x) + 1i*rand(sz,'like',grid.mesh.x)); % random initial guess
else
    phi = sqrt(nnn)*grid.normalize(phi0);
end

if(nargin<=4)
    tc=0;
end

ekk = exp(-grid.kk*dt);
% if(omega ~= 0)
%     ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt);
%     eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt);
% end    
MU = zeros(1000,1,'like',grid.mesh.x);
MU2 = zeros(1000,1,'like',grid.mesh.x);
EE = zeros(1000,1,'like',grid.mesh.x);
i = 0;
mgstep=10;

tmp2 = real(phi.*conj(phi)).*g+task.Fi(1:end-1,1:end-1,1:end-1)+V;
while true
    i=i+1;

%     if(omega ~= 0)
%         phi = grid.ifftx(ekx.*grid.fftx(phi));
%         phi = grid.iffty(eky.*grid.ffty(phi));
%     else
%         phi = grid.ifft(ekk.*grid.fft(phi));
%     end
    phi = exp(-tmp2*dt*0.5).*phi;
    phi = grid.ifft(ekk.*grid.fft(phi));
    if(omega ~= 0)
        lphi = phi;
        for ii = 1:3
            lphi = phi + dt*omega.*grid.lz(lphi);
            lphi = 0.5*(phi+lphi);
        end
        phi = phi + dt*omega.*grid.lz(lphi);
    end
    phi = exp(-tmp2*dt*0.5).*phi;
%     if(omega ~= 0)
%         phi = grid.iffty(eky.*grid.ffty(phi));
%         phi = grid.ifftx(ekx.*grid.fftx(phi));
%     else
%         phi = grid.ifft(ekk.*grid.fft(phi));
%     end

    tmp = real(phi.*conj(phi));
    mu = sqrt(task.Ntotal/grid.integrate(tmp));
    MU(i) = log(mu)/dt;
    phi=phi*mu;
    if(tc>0)
        phi = abs(phi).*exp(tc*1i*angs);
    end
    tmp = tmp*mu^2;
    if(mod(i,mgstep) == 5)
        task.set_pot_bc(phi);
        [task.Fi,~]=task.V_cycle(task.Fi,tmp+task.bar_dens,h,sz(1));
    end
    tmp2 = tmp.*g+task.Fi(1:end-1,1:end-1,1:end-1)+V;
    task.current_state = phi;
subplot(1,2,1)
imagesc(abs(phi(:,:,128)));
subplot(1,2,2)
imagesc(real(task.Fi(:,:,128)));
% subplot(2,2,3)
% plot(EE);
% subplot(2,2,4)
% plot(abs(abs(task.Fi(:,128,128))));
drawnow;
    if(i>50 && mod(i,10) == 0)
        if(nargout >= 3)
            MU2(i) = real(grid.inner(phi,task.applyham(phi)));
    %         MU2(i) = res;
        end
        if(nargout >= 4)
            EE(i) = task.get_energy(phi)/task.Ntotal;
%             EE(i)=max(abs(abs(phi(:))-abs(phi0(:))))/dt;
        end        
%         delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt;
        delta = max(abs(abs(phi(:))-abs(phi0(:))))/dt;
        if(delta < eps)
            if (dt<eps || dt<1e-4)
                break;
            else
                dt = dt/1.5;
                mgstep=mgstep*2;
                ekk = exp(-grid.kk*dt);
%                 if(omega ~= 0)
%                     ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt);
%                     eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt);                
%                 end
            end
        end
    end
    phi0=phi;
    if(i>=50000)
        warning('Convergence not reached');
        break;
    end
end

    task.current_mu = MU(end);
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

task.init_state = phi;
end
