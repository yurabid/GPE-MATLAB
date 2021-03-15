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
NNN=NN0;
gam = task.gamma;
% n_cn = task.n_crank;
n_rec = task.n_recalc;
% tau = task.decay_rate;
start = task.current_iter;

dt = ddt*1i/(1+1i*gam);
ekk = exp(-grid.etot*dt*0.5).*grid.mask;
% ekkp = exp(-grid.etot*dt*0.5).*grid.mask;
% ekkm = exp(grid.etot*dt*0.5).*grid.mask;
ekkphi = exp(-grid.rdphimat*dt*0.5);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end
variance = sqrt(task.T*gam*ddt);
sz = size(phi);
% mu = real(grid.inner(phi,task.applyham(phi)))./NN0;
mu = task.mu_init;
dt_outer = ddt*niter_inner;
% main BIG cycle starts here
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2)-mu;
%         phi = ekkp.*phi;
        % main SMALL cycle starts here
        for i=1:n_rec
            noise = variance*(randn(sz,'like',VV) + 1i*randn(sz,'like',VV));
            phi = ekk.*(phi+noise);
            phir = grid.ifftr(phi);
            phir = ekkphi.*phir;
            phirzp = grid.ifftzphi(phir);
            phirzp = exp(-(VV+g*abs(phirzp).^2)*dt).*phirzp;
            phir = grid.fftzphi(phirzp);
            phir = ekkphi.*phir;
            phi = grid.fftr(phir);
            phi = ekk.*(phi);            
            
            
%             phir = grid.sp2grid(phi);
%             phir = exp(-(VV-mu+g*abs(phir.^2))*dt).*phir;
%             phi = grid.grid2sp(phir) + variance*(randn(sz,'like',VV) + 1i*randn(sz,'like',VV));
%             phi = ekk.*phi;
        end
%         phi = ekkm.*phi;
        ncur = grid.integrate(real(phi.*conj(phi)));
    end
    err=task.ext_callback(phi,j,time2,mu,ncur);
    if(strcmp(err,'TERM'))
        break;
    end
end

