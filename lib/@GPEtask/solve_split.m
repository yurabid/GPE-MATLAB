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
tic;
grid = task.grid;
g = task.g;
if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = grid.integrate(abs(task.init_state).^2);
end
NNN=NN0;
gam = task.gamma;
n_rec = task.n_recalc;
if(mod(niter_inner,n_rec) ~= 0)
    n_rec = 1;
end
tau = task.decay_rate;
start = task.current_iter;

dt = ddt*1i/(1+1i*gam);
task.set_kinop(dt);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

task.nlinop = real(phi.*conj(phi));
muc = real(grid.inner(phi,task.applyham(phi)))./NN0;
if(task.mu_init > 0)
    mu = task.mu_init;
else
    mu = muc;
end
dt_outer = ddt*niter_inner;
% main BIG cycle starts here
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2)-mu;
        task.nlinop = exp(( - VV - g.*task.nlinop)*dt*0.5);
        phi = task.nlinop.*phi;
        % main SMALL cycle starts here
        for i=1:n_rec
            phi = task.ssft_kin_step(phi,dt);
            task.nlinop = exp(( - VV - g.*phi.*conj(phi))*dt);
            phi = task.nlinop.*phi;
        end
        task.nlinop = exp((VV + g.*phi.*conj(phi))*dt*0.5);
        phi = task.nlinop.*phi;
        task.nlinop = real(phi.*conj(phi));
        
        ncur = grid.integrate(task.nlinop);
        if(gam>0 && task.mu_init == 0)
            if(tau >0)
                NNN = NN0*exp(-time2/tau);
            end                
            phi = phi*sqrt(NNN/ncur);
            task.nlinop = task.nlinop*(NNN/ncur);
            muc = real(grid.inner(phi,task.applyham(phi,time2)))/NNN;
            mu = muc;
        else
            muc = real(grid.inner(phi,task.applyham(phi,time2)))/ncur;
        end
    end
    task.ext_callback(phi,j,time2,muc,ncur);
    
end

end
