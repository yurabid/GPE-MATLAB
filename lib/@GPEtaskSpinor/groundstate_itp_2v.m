function [phi, varargout] = groundstate_itp_2v(task,dt,eps,x0,phi)
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
ncomp = size(task.g,1);
task.ncomp = ncomp;
angs = -(atan2(grid.mesh.y,grid.mesh.x-x0) + atan2(grid.mesh.y,grid.mesh.x+x0))+pi;
VV = task.getVtotal(0);
V = VV{1};
coupl = task.coupling;
g = task.g;
if(nargin <= 4)
    phi = cellfun(@(x) grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V))./sqrt(ncomp), cell(1,ncomp), 'UniformOutput', false);
end
ekk = exp(-grid.kk*dt*0.5);
if(omega ~= 0)
    ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt);
    eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt);
end    
MU = zeros(5000,1,'like',V);
MU2 = zeros(5000,1,'like',V);
dts = zeros(5000,1,'like',V);
EE = zeros(5000,1,'like',V);
delta = 1;
i = 1;
iswitch = 50;

tmp2 = cell(1,ncomp);
cang = angle(coupl);
cosom = cosh(dt*abs(coupl));
sinomm = sinh(dt*abs(coupl)).*exp(-1i*cang);
sinomp = sinh(dt*abs(coupl)).*exp(1i*cang);

while true
    for j =1:ncomp
        tmp2{j} = VV{j};
        for k = 1:ncomp
            tmp2{j} = tmp2{j} + g(j,k)*abs(phi{k}).^2;%.*conj(phi{k}));
        end
    end
%     if(task.nlincpl~=0)
%         tmp3 = coupl + task.nlincpl*(abs(phi{1}).^2+abs(phi{2}).^2);
%         cosom = cosh(dt*tmp3);
%         sinomm = sinh(dt*tmp3);
%         sinomp = sinomm;
%     end
    if(omega ~= 0)
        phi{1} = grid.ifftx(ekx.*grid.fftx(phi{1}));
        phi{1} = grid.iffty(eky.*grid.ffty(phi{1}));
        phi{2} = grid.ifftx(ekx.*grid.fftx(phi{2}));
        phi{2} = grid.iffty(eky.*grid.ffty(phi{2}));        
    else
        phi{1} = grid.ifft(ekk.*grid.fft(phi{1}));
        phi{2} = grid.ifft(ekk.*grid.fft(phi{2}));
    end
    tmp2{1} = exp(-tmp2{1}*dt).*(cosom.*phi{1} - sinomp.*phi{2});
    tmp2{2} = exp(-tmp2{2}*dt).*(cosom.*phi{2} - sinomm.*phi{1});
    if(omega ~= 0)
        phi{1} = grid.ifftx(ekx.*grid.fftx(tmp2{1}));
        phi{1} = grid.iffty(eky.*grid.ffty(phi{1}));
        phi{2} = grid.ifftx(ekx.*grid.fftx(tmp2{2}));
        phi{2} = grid.iffty(eky.*grid.ffty(phi{2}));        
    else
        phi{1} = grid.ifft(ekk.*grid.fft(tmp2{1}));
        phi{2} = grid.ifft(ekk.*grid.fft(tmp2{2}));
    end
    
    ntot = real(grid.integrate(phi{1}.*conj(phi{1}) + phi{2}.*conj(phi{2})));
    mu = sqrt(task.Ntotal/ntot);
    MU(i) = log(mu)/dt;
    dts(i) = dt;
%     for j =1:ncomp
    phi{1}=abs(phi{1}).*exp(1i*(angs-angle(phi{2}))).*mu;
    phi{2}=phi{2}.*mu;
%     end
    task.current_state = phi;
    
    if(nargout >= 3)
        MU2(i) = real(task.inner(phi,task.applyham(phi)))/task.Ntotal;
    else
        MU2(i) = MU(i);
    end
    if(nargout >= 4)
        EE(i) = task.get_energy(phi)/task.Ntotal;
    else
        EE(i) = MU2(i);
    end

    if((i-iswitch)>200 && mod(i,10) == 0)
        delta = (abs(EE(i)-EE(i-9))/9 + abs(EE(i)-EE(i-1)))/dt;
        if(delta < eps)
            break;
        end
%         if(delta<dt*0.001 || delta < eps)
%             if (dt<eps || dt<5e-6)
%                 break;
%             else
%                 dt = dt/1.5;
%                 ekk = exp(-grid.kk*dt*0.5);
%                 if(omega ~= 0)
%                     ekx = exp(-(grid.kx.^2-2*grid.kx.*grid.mesh.y*omega)/4*dt);
%                     eky = exp(-(grid.ky.^2+2*grid.ky.*grid.mesh.x*omega)/4*dt);
%                 end                 
%                 cosom = cosh(dt*abs(coupl));
%                 sinomm = sinh(dt*abs(coupl)).*exp(-1i*cang);
%                 sinomp = sinh(dt*abs(coupl)).*exp(1i*cang);
%                 iswitch = i;
%             end
%         end
    end
    if(i>=500000)
        warning('Convergence not reached');
        break;
    end
    i=i+1;
end

if(nargout >= 2)
    MU = MU(1:i-1);
    varargout{1} = MU;
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
