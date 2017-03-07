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
ncomp = size(task.g,1);
coupl = task.coupling;
grid = task.grid;
% VV = task.getVtotal(0);
g = task.g;
gam = task.gamma;
n_rec = task.n_recalc;
tau = task.decay_rate;
start = task.current_iter;

if(start>0)
    phi = task.current_state;
else
    phi = task.init_state;
end

if(task.Ntotal > 0)
    NN0 = task.Ntotal;
else
    NN0 = 0;
    for j =1:ncomp
        NN0 = NN0 + grid.integrate(real(phi{j}.*conj(phi{j})));
    end
    NN0 = grid.integrate(abs(task.init_state).^2);
end
NNN=NN0;
ncur=NN0;

dt = ddt*1i/(1+1i*gam);
ekk = exp(-grid.kk*dt*0.5);

tmp = cell(1,ncomp);
tmp2 = cell(1,ncomp);
tmp3 = cell(1,ncomp);
cosom = cosh(dt*0.5*coupl);
sinom = sinh(dt*0.5*coupl);
hphi = task.applyham(phi);
mu = 0;
for j =1:ncomp
    mu = mu + real(grid.inner(phi{j},hphi{j}));
end
mu = mu/NN0;
dt_outer = ddt*niter_inner;
% main BIG cycle starts here
for j=start+1:niter_outer
    time=(j-1)*dt_outer;
    for jj=1:niter_inner/n_rec
        
        time2=time+(jj-1)*ddt*n_rec;
        VV = task.getVtotal(time2);

        % main SMALL cycle starts here
        for ii=1:n_rec

            for i =1:ncomp
                phi{i} = grid.ifft(ekk.*grid.fft(phi{i}));
            end
            for i =1:ncomp                
                tmp3{i} = cosom.*phi{i} - sinom.*phi{3-i};
            end
            for i =1:ncomp
                tmp2{i} = VV{i}-mu;
                for k = 1:ncomp
                    tmp2{i} = tmp2{i} + g(i,k)*real(tmp3{k}.*conj(tmp3{k}));
                end
            end             
            for i =1:ncomp                
                phi{i} = exp(-tmp2{i}*dt).*tmp3{i};
            end  
            for i =1:ncomp                
                tmp3{i} = cosom.*phi{i} - sinom.*phi{3-i};
            end            
            for i =1:ncomp
                phi{i} = grid.ifft(ekk.*grid.fft(tmp3{i}));
            end
               
        end
        
        if(gam>0)
            if(tau >0)
                NNN = NN0*exp(-time2/tau);
            end
            ncur = 0;
            mu = 0;
            for i =1:ncomp
                tmp{i} = real(phi{i}.*conj(phi{i}));
                ncur = ncur + grid.integrate(tmp{i});
            end  
            for i =1:ncomp
                phi{i} = phi{i}*sqrt(NNN/ncur);
            end
            hphi = task.applyham(phi);
            for i =1:ncomp
                mu = mu + real(grid.inner(phi{i},hphi{i}));
            end
            mu = mu/NNN;

        else
            ncur = 0;
            mu = 0;
            hphi = task.applyham(phi);
            for i =1:ncomp
                tmp{i} = real(phi{i}.*conj(phi{i}));
                ncur = ncur + grid.integrate(tmp{i});
                mu = mu + real(grid.inner(phi{i},hphi{i}));
            end
            mu = mu/NNN;
        end
    end
    task.ext_callback(phi,j,time2,mu,ncur);
    
end

end
