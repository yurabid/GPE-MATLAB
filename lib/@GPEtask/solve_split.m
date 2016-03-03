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

task.dispstat('','init');
% Initialize basic parameters
grid = task.grid;
% VV = task.getVtotal(0);
g = task.g;
omega = task.omega;
NN0=task.Ntotal;
NNN=NN0;
gam = task.gamma;
n_cn = task.n_crank;
n_rec = task.n_recalc;
tau = task.decay_rate;
start = task.current_iter;

dt = ddt*1i/(1+1i*gam);
ekk = exp(-grid.kk*dt);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

tmp2 = (phi.*conj(phi));

mu = real(grid.inner(phi,task.applyham(phi)))./NN0;

% muprev = mu;
% mu0=mu;
% dt = ddt*1i/(1+1i*gam); % time step (with gamma included)
dt_outer = ddt*niter_inner;
% ncur = grid.integrate(tmp2);
% main BIG cycle starts here
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        
        mu_run = mu;
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2);
        phi = exp((mu - VV - g*tmp2)*dt/2).*phi;
        % main SMALL cycle starts here
        for i=1:n_rec
            phi = grid.ifft(ekk.*grid.fft(phi));
            if(omega ~= 0)
                lphi = phi;
                for ii = 1:n_cn
                    lphi = phi + dt*omega*grid.lz(lphi);
                    lphi = 0.5*(phi+lphi);
                end
                phi = phi + dt*omega*grid.lz(lphi);
            end
            
            %         if(tau>0)
            %             mu_run = (mu + (mu0-mu)*i*ddt/dt_outer)*exp(-i*ddt/tau);
            %             mu_run = mu*exp(-i*ddt/tau);
            %         end
            phi = exp((mu_run - VV - g*phi.*conj(phi))*dt).*phi;
        end
        
        phi = exp((VV - mu + g*phi.*conj(phi))*dt/2).*phi;
        tmp2 = (phi.*conj(phi));
        
        if(gam>0)
            if(tau >0)
                NNN = NN0*exp(-time2/tau);
            end
            
            %         nprev = ncur;
            ncur = grid.integrate(tmp2);
            phi = phi*sqrt(NNN/ncur);
            tmp2 = tmp2*(NNN/ncur);
            mu = real(grid.inner(phi,task.applyham(phi)))/NNN;
            %         muprev = (mu0+muprev)/2;
            %         mu = mu - log(ncur/NNN)/(dt_outer*2*gam);
            %         mu = real(mu0) - log(ncur/NNN)/(dt_outer*2*gam);% - (ncur-nprev)/(gam*(ncur+nprev)*dt_outer);

        end
    end
    task.ext_callback(phi,j,time2,mu,NNN);
    
end

end
