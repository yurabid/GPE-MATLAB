function [phi, varargout] = groundstate_besp(task,dt,eps,phi)
% groundstate_itp - Calculate the stationary state of 2-component GPE 
% with BESP Imaginary Time Propagation method.

% Based on:
% Bao, W., Chern, I.L. and Lim, F.Y., 2006. Efficient and spectrally accurate 
% numerical methods for computing ground and first excited states in 
% Bose–Einstein condensates. Journal of Computational Physics, 219(2), pp.836-854.


%  Input
%    dt    :  evolution time step
%    eps   :  desired accuracy (applied to chemical potential)
%    phi0  :  (optional) initial approximation of the wave function
%  Output
%    phi      :  calculated stationary state
%    mu       :  (optional) history of chemical potential values during evolution
%    mu2      :  (optional) history of chemical potential from integral evaluation
%    E        :  (optional) history of total energy values

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
MU = zeros(5000,1,'like',V);
MU2 = zeros(5000,1,'like',V);
EE = zeros(5000,1,'like',V);
delta = 1;
i = 1;
iswitch = 20;
eps2=eps*100;
tmp2 = cell(1,ncomp);

while true
    tmp2{1} = VV{1} + g(1,1)*abs(phi{1}).^2 + g(1,2)*abs(phi{2}).^2;%.*conj(phi{k}));
    tmp2{2} = VV{2} + g(2,2)*abs(phi{2}).^2 + g(2,1)*abs(phi{1}).^2;%.*conj(phi{k}));

	bmax1 = max(tmp2{1}(:));
	bmin1 = min(tmp2{1}(:));
    bmax2 = max(tmp2{2}(:));
	bmin2 = min(tmp2{2}(:));    
	alpha = (bmax1+bmin1+bmax2+bmin2)/4;
	phi1hat = grid.fft(phi{1});
    phi2hat = grid.fft(phi{2});
    
	g1hat = grid.fft((alpha-tmp2{1}).*phi{1}-coupl.*phi{2});
    g2hat = grid.fft((alpha-tmp2{2}).*phi{2}-conj(coupl).*phi{1});
	phi1 = grid.ifft((phi1hat+dt*g1hat)./(1+dt*(alpha+grid.kk)));
    phi2 = grid.ifft((phi2hat+dt*g2hat)./(1+dt*(alpha+grid.kk)));
    
    delta=1;
    while delta>eps2
        g1hat = grid.fft((alpha-tmp2{1}).*phi1-coupl.*phi2);
        g2hat = grid.fft((alpha-tmp2{2}).*phi2-conj(coupl).*phi1);
        phi11 = grid.ifft((phi1hat+dt*g1hat)./(1+dt*(alpha+grid.kk)));
        phi22 = grid.ifft((phi2hat+dt*g2hat)./(1+dt*(alpha+grid.kk)));
        delta = max(abs([phi1(:)-phi11(:);phi2(:)-phi22(:)]));
        phi1 = phi11;
        phi2 = phi22;
    end

    ntot = real(grid.integrate(phi1.*conj(phi1) + phi2.*conj(phi2)));
    mu = sqrt(task.Ntotal/ntot);
    phi1=phi1.*mu;
    phi2=phi2.*mu;
    
    MU(i) = (mu-1)/dt;  
    
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
    if i==iswitch
        eps2=eps;
    end
    if((i-iswitch)>0 && mod(i,5) == 0)
%         delta = (abs(EE(i)-EE(i-9))/9 + abs(EE(i)-EE(i-1)))/dt;
%         delta = abs((EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)))/EE(i);
        delta = max(abs([phi{1}(:)-phi1(:);phi{2}(:)-phi2(:)]));
        if(delta < eps)
            break;
        end

    end
    phi{1} = phi1;
    phi{2} = phi2;
    if(i>=50000)
        warning('Convergence not reached');
        break;
    end
    i=i+1;
%     imagesc(abs(phi{1}(:,:,end/2)));drawnow;
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
