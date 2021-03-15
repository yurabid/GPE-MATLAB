function [phi, varargout] = bdg_itp(task,dt,eps,num)
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
omega = task.omega;
phi0 = task.current_state;
NN = grid.inner(phi0,phi0);
ekk = exp(-grid.kk*0.5*dt);
if(omega ~= 0)
    ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*task.omega)/4*dt);
    eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*task.omega)/4*dt);
end
MU = zeros(1000,1,'like',V);
% MU2 = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);
task.spectrum_u = zeros([size(phi0),num]);
task.spectrum_v = zeros([size(phi0),num]);
task.spectrum_w = zeros(num,1);
mu = real(grid.inner(phi,task.applyham(phi)));
tmp = real(phi0.*conj(phi0)).*g;
tmp2 = 2*tmp+V;
cosnl = cosh(dt*tmp);
sinnl = sinh(dt*tmp);

for j=1:num
uc = sqrt(nnn)*grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V));
vc = sqrt(nnn)*grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V));

i = 0;
while true
    i=i+1;
    if(omega ~= 0)
        uc = grid.ifftx(ekx.*grid.fftx(uc));
        uc = grid.iffty(eky.*grid.ffty(uc));
        vc = grid.ifftx(ekx.*grid.fftx(vc));
        vc = grid.iffty(eky.*grid.ffty(vc));
    else
        uc = grid.ifft(ekk.*grid.fft(uc));
        vc = grid.ifft(ekk.*grid.fft(vc));
    end
    uc = exp(-tmp2*dt).*(cosnl.*uc-sinnl.*vc);
    vc = exp(-tmp2*dt).*(cosnl.*vc-sinnl.*uc);
    if(omega ~= 0)
        uc = grid.ifftx(ekx.*grid.fftx(uc));
        uc = grid.iffty(eky.*grid.ffty(uc));
        vc = grid.ifftx(ekx.*grid.fftx(vc));
        vc = grid.iffty(eky.*grid.ffty(vc));
    else
        uc = grid.ifft(ekk.*grid.fft(uc));
        vc = grid.ifft(ekk.*grid.fft(vc));
    end

tmp3 = grid.inner(vc,phi0)*phi0/NN;
tmp4 = grid.inner(uc,phi0)*phi0/NN;
if (j>1)
    for jj=1:j-1
        tmp3 = tmp3 + grid.inner(vc,task.spectrum_v(:,:,:,jj))*task.spectrum_v(:,:,:,jj);
        tmp4 = tmp4 + grid.inner(uc,task.spectrum_u(:,:,:,jj))*task.spectrum_u(:,:,:,jj);
    end
end
    vc = vc - tmp3;
    uc = uc - tmp4;
    mu = sqrt(1/(grid.integrate(abs(vc).^2+abs(uc).^2)));
    MU(i) = log(mu)/dt;
    vc=vc*mu;
    uc=uc*mu;
%     tmp = tmp*mu^2;
%     tmp2 = tmp.*g+V;
%     task.current_state = phi;
%     imagesc(abs(phi));drawnow;
%     if(nargout >= 3)
%         MU2(i) = real(grid.inner(phi,task.applyham(phi)));
%     end
%     if(nargout >= 4)
%         EE(i) = task.get_energy(phi)/task.Ntotal;
%     end

    if(i>50 && mod(i,10) == 0)
        delta = (abs(MU(i)-MU(i-9))/9 + abs(MU(i)-MU(i-1)))/dt/MU(i);
        if(delta < eps)
            if (dt<eps*10 || dt<1e-4)
                break;
            else
                dt = dt/1.5;
                ekk = exp(-grid.kk*0.5*dt);
                if(omega ~= 0)
                    ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*task.omega)/4*dt);
                    eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*task.omega)/4*dt);
                end
            end
        end
    end

    if(i>=50000)
        warning('Convergence not reached');
        break;
    end
end

if(nargout >= 2)
    MU = MU(1:i);
    if(task.Ntotal > 0)
        MUEX = MU(i) - (MU(i)-MU(i-5))^2/(MU(i)-2*MU(i-5)+MU(i-10)); % exponential extrapolation
        MU = [MU; MUEX];
        task.current_mu = MUEX;
        task.current_n = task.Ntotal;
    end
    varargout{1} = MU;
end
if(nargout >= 3)
    if(task.Ntotal > 0)
        MUEX = MU2(i) - (MU2(i)-MU2(i-5))^2/(MU2(i)-2*MU2(i-5)+MU2(i-10)); % exponential extrapolation
        MU2 = [MU2(1:i); MUEX]/task.Ntotal;
    else
        MU2 = MU2(1:i)./MU;
        task.current_mu = MU2(end);
        task.current_n = MU(end);
    end
    varargout{2} = MU2;
end
if(nargout >= 4)
    if(task.Ntotal > 0)
        MUEX = EE(i) - (EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)); % exponential extrapolation
        EE = [EE(1:i); MUEX];
    else
        EE = EE(1:i)./MU;
    end
    varargout{3} = EE;
end

task.init_state = phi;
end
