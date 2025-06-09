function phi = solve_split(task,ddt,niter_inner,niter_outer)
% Calculate the dynamics of GPE-Poisson system 
% using split-step Fourier method. The gravitational potential is
% recalculated at each time step with a multigrid V-cycle
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
sz = grid.nx+1;
qms = zeros(9,5,'like',grid.x);
start = task.current_iter;
dt = ddt*1i/(1+1i*task.gamma);
task.set_kinop(dt);
if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end
if(numel(task.Fi)==0)
    task.set_pot(abs(phi).^2,true);
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
        task.nlinop = exp(-0.5*dt*(task.g.*tmp+task.getVtotal(time2)-mu));
        phi1 = task.nlinop.*phi;
        phi1 = task.ssft_kin_step(phi1,dt);
        phi1 = task.nlinop.*phi1;

        tmp=(tmp+real(phi1.*conj(phi1)))/2;
        [task.Fi,~]=task.V_cycle(task.Fi,tmp,grid.dx,sz); 
        
        task.nlinop = exp(-0.5*dt*(task.g.*tmp+task.getVtotal(time2)-mu));
        phi = task.nlinop.*phi;
        phi = task.ssft_kin_step(phi,dt);
        phi = task.nlinop.*phi;
        
        tmp = real(phi.*conj(phi));
        task.set_pot_bc(tmp);
        % [task.Fi,~]=task.V_cycle(task.Fi,tmp,h,sz);

        if(task.calculate_qmder)
            qms = [qms(:,2:5),task.current_qm(:)];
            task.current_qmder = qms*[-0.5;1;0;-1;0.5]/ddt^3;
        end

    end
    ncur = grid.integrate(tmp);
    muc = real(grid.inner(phi,task.applyham(phi,time2)))/ncur;
    task.ext_callback(phi,j,time2,muc,ncur);    
end

end
