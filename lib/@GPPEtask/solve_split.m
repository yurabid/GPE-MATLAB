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
% VV = task.getVtotal(0);
g = task.g;
h=grid.x(2)-grid.x(1);
sz = size(grid.mesh.x)+1;
if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = grid.integrate(abs(task.init_state).^2);
end
gam = task.gamma;
qms = zeros(9,5,'like',grid.x);
start = task.current_iter;

dt = ddt*1i/(1+1i*gam);
ekk = exp(-grid.kk*dt);

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

tmp2 = real(phi.*conj(phi));
muc = real(grid.inner(phi,task.applyham(phi)))./NN0;
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
        VV = g.*tmp2+task.getVtotal(time2)-mu+task.Fi(1:end-1,1:end-1,1:end-1);
        phi1 = exp(-VV*dt*0.5).*phi;
        phi1 = grid.ifft(ekk.*grid.fft(phi1));
        phi1 = exp(-VV*dt*0.5).*phi1;

        tmp2=(tmp2+real(phi1.*conj(phi1)))/2;
        [task.Fi,~]=task.V_cycle(task.Fi,tmp2+task.bar_dens,h,sz(1)); 
        
        VV = g.*tmp2+task.getVtotal(time2)-mu+task.Fi(1:end-1,1:end-1,1:end-1);
        phi = exp(-VV*dt*0.5).*phi;
        phi = grid.ifft(ekk.*grid.fft(phi));
        phi = exp(-VV*dt*0.5).*phi;
        
        tmp2 = real(phi.*conj(phi));
        task.set_pot_bc(phi);
        [task.Fi,~]=task.V_cycle(task.Fi,tmp2+task.bar_dens,h,sz(1));
        qms = [qms(:,2:5),task.current_qm(:)];
        task.current_qmder = qms*[-0.5;1;0;-1;0.5]/ddt^3;
        
        ncur = grid.integrate(tmp2);
        muc = real(grid.inner(phi,task.applyham(phi,time2)))/ncur;

    end
    task.ext_callback(phi,j,time2,muc,ncur);    
end

end
