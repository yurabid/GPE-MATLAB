function phi = solve_split_ft(task,ddt,niter_inner,niter_outer)
% solve_split - Calculate the dynamics of GPPE with the split-step method.
% Gravitational potential is calculated with FFT method
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
kk_inv = 1./grid.kk/2;
kk_inv(1,1,1) = 0;
g = task.g;
start = task.current_iter;
dt = ddt*1i/(1+1i*task.gamma);
task.set_kinop(dt);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

tmp = real(phi.*conj(phi));
muc = real(grid.inner(phi,task.applyham(phi)))./grid.integrate(tmp);
if(task.mu_init > 0)
    mu = task.mu_init;
else
    mu = muc;
end
dt_outer = ddt*niter_inner;
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner
        time2=time+(jj-1)*ddt;
        task.Fi = -grid.ifft(kk_inv.*grid.fft(task.nlinop));
        task.nlinop = exp(-0.5*dt*(g.*task.nlinop+task.getVtotal(time2)-mu));
        phi = task.nlinop.*phi;
        phi = task.ssft_kin_step(phi,dt);
        phi = task.nlinop.*phi;
        task.nlinop = real(phi.*conj(phi));
    end
    ncur = grid.integrate(task.nlinop);
    muc = real(grid.inner(phi,task.applyham(phi,time2)))/ncur;
    task.ext_callback(phi,j,time2,muc,ncur);    
end

end
