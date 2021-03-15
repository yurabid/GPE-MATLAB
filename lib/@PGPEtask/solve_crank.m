function phi = solve_crank(task,ddt,niter_inner,niter_outer)
% solve_split - Calculate the dynamics of GPE with the semi-implicit Ckank-Nicholson scheme.
%
%  Usage :
%    phi = task.solve_crank(ddt,niter_inner,niter_outer)
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

phir = grid.sp2grid(phi);
variance = sqrt(task.T*gam*ddt);
sz = size(phi);
% mu = real(grid.inner(phi,task.applyham(phi,0,phir)))./NN0;
mu = task.mu_init;
dt_outer = ddt*niter_inner;
% main BIG cycle starts here
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2);

        % main SMALL cycle starts here
        for i=1:n_rec

            lphi = phi;
            lphir = phir;
            noise = variance*(randn(sz,'like',VV) + 1i*randn(sz,'like',VV));
            for ii = 1:n_cn        
                lphi = phi - dt*((grid.etot-mu).*lphi + grid.grid2sp((VV+g*abs(lphir.^2)).*lphir));
                lphi = 0.5*(phi+lphi);
                lphir = grid.sp2grid(lphi);
            end
            phi = (phi - dt*((grid.etot-mu).*lphi + grid.grid2sp((VV+g*abs(lphir.^2)).*lphir)) - noise).*grid.mask;
            phir = grid.sp2grid(phi);
        end
        
        tmp2 = real(phi.*conj(phi));
        
%         if(gam>0)
%             if(tau >0)
%                 NNN = NN0*exp(-time2/tau);
%             end
%             
%             ncur = grid.integrate(tmp2);
%             phi = phi*sqrt(NNN/ncur);
%             phir = phir*sqrt(NNN/ncur);
%             mu = real(grid.inner(phi,task.applyham(phi,time2,phir)))/NNN;
% 
%         else
            ncur = grid.integrate(tmp2);
%             mu = real(grid.inner(phi,task.applyham(phi,time2,phir)))/NNN;
%         end
    end
    task.ext_callback(phi,j,time2,mu,ncur);
end

end
