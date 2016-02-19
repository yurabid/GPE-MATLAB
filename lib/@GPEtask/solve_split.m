function phi = solve_split(task,ddt,niter_inner,niter_outer)
% solve_split - Calculate the dynamics of GPE with the split-step method.
%
%  Usage :
%    phi = task.solve_split(ddt,niter_inner,niter_outer)
%  Input
%    ddt         :  evolution time step
%    niter_inner :  number of internal iterations between callbacks
%    niter_outer :  number of external iterations
%  Output
%    phi       :  final state


% Initialize basic parameters
grid = task.grid;
VV = task.getVtotal(0);
g = task.g;
omega = task.omega;
NN0=task.Ntotal;
gam = task.gamma;
n_cn = 10;
tau = task.decay_rate;
start = task.current_iter;

dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
ekk = exp(-grid.kk*dt);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

tmp2 = abs(phi.*conj(phi));

% mu = abs(grid.integrate(abs(conj(phi).*grid.ifft(grid.kk.*grid.fft(phi))) + (VV+g*tmp2).*tmp2))./NN0 - omega*LLL;
mu = grid.inner(phi,task.applyham(phi))./NN0;
if(omega ~= 0)
    LLL = -imag(grid.integrate(conj(phi).*(grid.mesh.x.*grid.derivy(phi) -...
        grid.mesh.y.*grid.derivx(phi))))./NN0;
    mu = mu -omega*LLL;
end
muprev = mu;
dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
dt_outer = ddt*niter_inner;

% main BIG cycle starts here
for j=start+1:niter_outer

    time=(j-1)*dt_outer;
    phi = exp(-(VV - mu  + g*tmp2)*dt/2).*phi;
    mu_run = mu;
    % main SMALL cycle starts here
    for i=1:niter_inner
        phi = grid.ifft(ekk.*grid.fft(phi));
        if(omega ~= 0)
            lphi = phi;
            for ii = 1:n_cn
                lphi = phi + dt*1i*omega*(grid.mesh.x.*grid.derivy(lphi) - ...
                    grid.mesh.y.*grid.derivx(lphi));
                lphi = 0.5*(phi+lphi);
            end
            phi = phi + dt*1i*omega*(grid.mesh.x.*grid.derivy(lphi) - ...
                grid.mesh.y.*grid.derivx(lphi));
        end
        time2=time+(i-1)*ddt;
        VV = task.getVtotal(time2);
        if(tau>0)
            mu_run = mu*exp(-i*ddt/tau);
        end
        phi = exp((mu_run - VV - g*phi.*conj(phi))*dt).*phi;
    end
    
    phi = exp((VV - mu + g*phi.*conj(phi))*dt/2).*phi;
    tmp2 = abs(phi.*conj(phi));

    if(gam>0 && tau >0)
        NNN = NN0*exp(-time2/tau);
        muprev = (mu+muprev)/2;
        mu = muprev + (-log(grid.integrate(tmp2)/NNN)/(dt_outer*2*gam));
    end
    
    task.ext_callback(phi,j,time2,mu);
    
end

end
