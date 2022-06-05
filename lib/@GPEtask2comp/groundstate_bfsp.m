function [phi, varargout] = groundstate_bfsp(task,dt,eps,phi)
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
MU = zeros(5000,1,'like',V);
MU2 = zeros(5000,1,'like',V);
dts = zeros(5000,1,'like',V);
EE = zeros(5000,1,'like',V);
i = 1;
iswitch = 50;

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
    
    ntot = real(grid.integrate(phi1.*conj(phi1) + phi2.*conj(phi2)));
    mu = sqrt(task.Ntotal/ntot);
    phi1=phi1.*mu;
    phi2=phi2.*mu;
    

    MU(i) = log(mu)/dt;
    dts(i) = dt;    
    
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

    if((i-iswitch)>50 && mod(i,10) == 0)
%         delta = (abs(EE(i)-EE(i-9))/9 + abs(EE(i)-EE(i-1)))/dt;
%         delta = abs((EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)))/EE(i);
        delta = max(abs([phi{1}(:)-phi1(:);phi{2}(:)-phi2(:)]));
        if(delta < eps)
            if (dt<eps*10)
                break;
            else
                dt = dt/3; 
                iswitch = i;
            end
        end

    end
%     imagesc(abs(phi1));drawnow;
%     phi1 = [phi1(1:end,1:end/2),conj(flip(phi1(1:end,1:end/2),2))];
%     phi1 = [phi1(1:end/2,1:end/2),conj(flip(phi1(1:end/2,1:end/2),2));conj(flip(phi1(1:end/2,1:end/2),1)),flip(flip(phi1(1:end/2,1:end/2),1),2)];
    phi{1} = phi1;
    phi{2} = phi2;
    if(i>=50000)
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
