function txt=post_process(task)
    global vv tshift rscale;
%     curx = vv*task.current_time;
    xm = 0.7e-6/rscale;
    curx = xm-bar_pos(task.current_time,xm,tshift);
    phi = task.current_state;
    dx = task.grid.x(2)-task.grid.x(1);
    w2d = dx*dx;
    ndim = numel(size(task.grid.mesh.x));
    ddens = 2*imag(conj(phi).*(task.applyham(phi)-task.current_mu*phi)/(1+1i*task.gamma)); % - phi.*conj(task.applyham(phi)));
    [z2,x2] = meshgrid(task.grid.z,task.grid.x);
    if(ndim==3)
        slice = squeeze(phi(end/2,:,:));
        densz = sum(abs(phi).^2,3)*(task.grid.z(2)-task.grid.z(1));
        ddensz = sum((ddens),3)*(task.grid.z(2)-task.grid.z(1));
        mulocal = sum(conj(phi).*task.applyham(phi),3)*(task.grid.z(2)-task.grid.z(1));
    else
        slice = phi;
        densz = abs(phi).^2;
        ddensz = (ddens);
        mulocal = conj(phi).*task.applyham(phi);
    end

%     LLL = -imag(task.grid.integrate(conj(phi).*(task.grid.mesh.x.*grid.derivy(phi) -...
%         task.grid.mesh.y.*task.grid.derivx(phi))))./task.current_n;

    MU = real(sum(sum(mulocal))*w2d/task.current_n);
    
%     ang1 = angle(sum(sum(slice.*(x2<curx))));
%     dphase = angle(sum(sum(slice.*exp(-1i*ang1).*(x2>curx))));
    ang1 = angle(sum(sum(slice(1:end/2,:))));
    dphase = angle(sum(sum(slice(end/2+1:end,:).*exp(-1i*ang1)))); 
    ang1 = angle(sum(slice(1:end/2,:),1));
    ang2 = angle(sum(slice(end/2+1:end,:),1));
    [coresp, coresm] = detect_core(slice,z2,x2);
    N1 = sum(sum(densz.*(task.grid.mesh.x2<curx)))*w2d;
    N2 = sum(sum(densz.*(task.grid.mesh.x2>curx)))*w2d;  
    
    dN1 = sum(sum(ddensz.*(task.grid.mesh.x2<curx)))*w2d;
    dN2 = sum(sum(ddensz.*(task.grid.mesh.x2>curx)))*w2d;

    mu1 = sum(sum(mulocal.*(task.grid.mesh.x2<curx)))*w2d/N1;
    mu2 = sum(sum(mulocal.*(task.grid.mesh.x2>curx)))*w2d/N2;
    
    dmu = (mu1-mu2); 
    task.history.N1(task.current_iter) = N1;
    task.history.N2(task.current_iter) = N2;
    task.history.dN1(task.current_iter) = dN1;
    task.history.dN2(task.current_iter) = dN2;
    task.history.dmu(task.current_iter) = dmu;
    task.history.dphase(task.current_iter) = dphase;
    task.history.mu_int(task.current_iter) = MU;
%     plot(task.grid.z,(ang1-ang2)); hold on; 
%     plot(task.grid.z,log(abs(tan((ang1-ang2)/4)))); ylim([-7,7]);
%     save(sprintf('snapshots/slice_%05d',task.current_iter),'slice','mulocal','densz','curx','N1','N2','dmu','dphase','MU');
    imagesc(task.grid.z,task.grid.x,squeeze(abs(phi(end/2,:,:)))); hold on;
    plot(coresm(:,1),coresm(:,2),'MarkerSize',4,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1);
    plot(coresp(:,1),coresp(:,2),'MarkerSize',4,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1);
    drawnow; hold off;
    txt = sprintf('z = %0.3f', (N1-N2)/task.Ntotal);
end
