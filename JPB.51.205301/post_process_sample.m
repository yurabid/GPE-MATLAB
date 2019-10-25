function txt=post_process_sample(task)
    global vv;
    curx = vv*task.current_time; % current position of the barrier
    phi = task.current_state;
    w2d = (task.grid.x(2)-task.grid.x(1))*(task.grid.y(2)-task.grid.y(1));
    slice = squeeze(phi(end/2,:,:));
    densz = sum(abs(phi).^2,3)*(task.grid.z(2)-task.grid.z(1));
    mulocal = sum(conj(phi).*task.applyham(phi),3)*(task.grid.z(2)-task.grid.z(1));
    MU = real(sum(sum(mulocal))*w2d/task.current_n);
    
    ang1 = angle(sum(sum(slice.*(task.grid.mesh.x2<curx))));
    dphase = angle(sum(sum(slice.*exp(-1i*ang1).*(task.grid.mesh.x2>curx))));
    N1 = sum(sum(densz.*(task.grid.mesh.x2<curx)))*w2d;
    N2 = sum(sum(densz.*(task.grid.mesh.x2>curx)))*w2d;  
    mu1 = sum(sum(mulocal.*(task.grid.mesh.x2<curx)))*w2d/N1;
    mu2 = sum(sum(mulocal.*(task.grid.mesh.x2>curx)))*w2d/N2;
    dmu = (mu1-mu2); 
    task.history.Z(task.current_iter) = (N1-N2)/task.Ntotal; % population imballance
    task.history.dmu(task.current_iter) = dmu; % chemical potential difference
    task.history.dphase(task.current_iter) = dphase; % phase difference
    task.history.mu_int(task.current_iter) = MU; % total chemical potential
    txt = sprintf('z = %0.3f', (N1-N2)/task.Ntotal);
end
