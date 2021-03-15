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
grid = task.grid;
g = task.g;
omega = task.omega;
if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = grid.integrate(abs(task.init_state).^2);
end
% NNN=NN0;
gam = task.gamma;
n_cn = task.n_crank;
n_rec = task.n_recalc;
tau = task.decay_rate;
start = task.current_iter;

dt = ddt*1i/(1+1i*gam);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

tmp2 = real(phi.*conj(phi));
% task.mu_init = real(grid.inner(task.init_state,task.applyham(task.init_state)))./NN0;
if(task.ecut>0 && task.T>0)
    mask = grid.kk < task.ecut;
    ekk = exp(-grid.kk*dt).*mask; % mask in fourier space implements cut-off of high-energy modes 
else
    ekk = exp(-grid.kk*dt);
end
dt_outer = ddt*niter_inner;
sz = size(phi);
variance = sqrt(task.T*gam*ddt/grid.weight);
% main BIG cycle starts here
for j=start+1:niter_outer
    mu = task.mu_init;
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2);
        phi = exp((mu - VV - g*tmp2)*dt/2).*phi;
        % main SMALL cycle starts here
        for i=1:n_rec
%             phis = grid.fft(phi) + sqrt(grid.nx*grid.ny*grid.nz)*variance*(randn(sz,'like',grid.mesh.x) + 1i*randn(sz,'like',grid.mesh.x)).*mask;
            phi = phi + variance*(randn(sz,'like',VV) + 1i*randn(sz,'like',VV));
            phi = grid.ifft(ekk.*grid.fft(phi));
            if(omega ~= 0)
                lphi = phi;
                for ii = 1:n_cn
                    lphi = phi + dt*omega*grid.lz(lphi);
                    lphi = 0.5*(phi+lphi);
                end
                phi = phi + dt*omega*grid.lz(lphi);
            end
            phi = exp((mu - VV - g*phi.*conj(phi))*dt).*phi;
        end
        
        phi = exp((VV - mu + g*phi.*conj(phi))*dt/2).*phi;
        tmp2 = real(phi.*conj(phi));
               
        ncur = grid.integrate(tmp2);
    end
    err=task.ext_callback(phi,j,time2,mu,ncur);
    if(strcmp(err,'TERM'))
        break;
    end
end
