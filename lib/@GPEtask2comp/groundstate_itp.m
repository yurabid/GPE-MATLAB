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
omega=task.omega;
% ncomp = size(task.g,1);
% task.ncomp = ncomp;
VV = task.getVtotal(0);
V = VV{1};
% coupl = task.coupling;
g = task.g;
if(nargin <= 3)
    phi = cellfun(@(x) grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V))./sqrt(2), cell(1,2), 'UniformOutput', false);
end
ekk1 = exp(-grid.kk*dt*0.5/task.M1);
ekk2 = exp(-grid.kk*dt*0.5/task.M2);
if(omega ~= 0)
    ekx1 = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt/task.M1);
    eky1 = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt/task.M1);
    ekx2 = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt/task.M2);
    eky2 = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt/task.M2);    
end    
MU1 = zeros(5000,1,'like',V);
MU2 = zeros(5000,1,'like',V);
dts = zeros(5000,1,'like',V);
EE = zeros(5000,1,'like',V);
% delta = 1;
i = 1;
iswitch = 20;

tmp2 = cell(1,2);
% cang = angle(coupl);
% cosom = cosh(dt*abs(coupl));
% sinomm = sinh(dt*abs(coupl)).*exp(-1i*cang);
% sinomp = sinh(dt*abs(coupl)).*exp(1i*cang);

while true
%     for j =1:ncomp
%         tmp2{j} = VV{j};
%         for k = 1:ncomp
%             tmp2{j} = tmp2{j} + g(j,k)*abs(phi{k}).^2;%.*conj(phi{k}));
%         end
%     end
    tmp2{1} = VV{1} + g(1,1)*abs(phi{1}).^2 + g(1,2)*abs(phi{2}).^2;
    tmp2{2} = VV{2} + g(2,1)*abs(phi{1}).^2 + g(2,2)*abs(phi{2}).^2;
    
    if(omega ~= 0)
        phi{1} = grid.ifftx(ekx1.*grid.fftx(phi{1}));
        phi{1} = grid.iffty(eky1.*grid.ffty(phi{1}));
        phi{2} = grid.ifftx(ekx2.*grid.fftx(phi{2}));
        phi{2} = grid.iffty(eky2.*grid.ffty(phi{2}));        
    else
        phi{1} = grid.ifft(ekk1.*grid.fft(phi{1}));
        phi{2} = grid.ifft(ekk2.*grid.fft(phi{2}));
    end
    tmp2{1} = exp(-tmp2{1}*dt).*phi{1} ;
    tmp2{2} = exp(-tmp2{2}*dt).*phi{2} ;
    if(omega ~= 0)
        phi{1} = grid.iffty(eky1.*grid.ffty(phi{1}));
        phi{1} = grid.ifftx(ekx1.*grid.fftx(tmp2{1}));
        phi{2} = grid.iffty(eky2.*grid.ffty(phi{2}));        
        phi{2} = grid.ifftx(ekx2.*grid.fftx(tmp2{2}));
    else
        phi{1} = grid.ifft(ekk1.*grid.fft(tmp2{1}));
        phi{2} = grid.ifft(ekk2.*grid.fft(tmp2{2}));
    end
    
    n1c = real(grid.integrate(phi{1}.*conj(phi{1})));
    n2c = real(grid.integrate(phi{2}.*conj(phi{2})));
    mu1 = sqrt(task.N1/n1c);
    mu2 = sqrt(task.N2/n2c);
    MU1(i) = log(mu1)/dt;
    MU2(i) = log(mu2)/dt;
    dts(i) = dt;
    phi{1}=phi{1}.*mu1;
    phi{2}=phi{2}.*mu2;
    task.current_state = phi;
    
%     if(nargout >= 3)
%         MU2(i) = real(task.inner(phi,task.applyham(phi)))/task.Ntotal;
%     else
%         MU2(i) = MU1(i);
%     end
    if(nargout >= 4)
        EE(i) = task.get_energy(phi)/task.Ntotal;
    else
        EE(i) = MU2(i);
    end
%     subplot(2,2,1);
%     imagesc(abs(phi{1}));
%     subplot(2,2,2);
%     imagesc(angle(phi{1}));
%     subplot(2,2,3);
%     imagesc(abs(phi{2}));
%     subplot(2,2,4);
%     plot([MU1(max(i-200,1):i),MU2(max(i-200,1):i)]);
%     imagesc(angle(phi{2}));
%     drawnow;
    if((i-iswitch)>10 && mod(i,10) == 0)
        delta = (abs(EE(i)-EE(i-9))/9 + abs(EE(i)-EE(i-1)))/dt;
%         if(delta < eps)
%             break;
%         end
        if(delta < eps)
            if (dt<eps || dt<5e-6)
                break;
            else
                dt = dt/1.5;
                ekk1 = exp(-grid.kk*dt*0.5/task.M1);
                ekk2 = exp(-grid.kk*dt*0.5/task.M2);
                if(omega ~= 0)
                    ekx1 = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt/task.M1);
                    eky1 = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt/task.M1);
                    ekx2 = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt/task.M2);
                    eky2 = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt/task.M2);    
                end                

                iswitch = i;
            end
        end
    end
    if(i>=500000)
        warning('Convergence not reached');
        break;
    end
    i=i+1;
end

if(nargout >= 2)
    MU1 = MU1(1:i-1);
    varargout{1} = MU1;
end
if(nargout >= 3)
    MU2 = MU2(1:i-1);
    varargout{2} = MU2;
end
if(nargout >= 4)
    EE = EE(1:i-1);
    varargout{3} = EE;
end
task.init_state = phi;
task.current_state = phi;
end
