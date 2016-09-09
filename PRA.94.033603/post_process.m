function txt=post_process(task)
    global omegaBar;
    angles3 = atan2(task.grid.mesh.x2,task.grid.mesh.y2);
    curangle = rotangle(task.current_time,omegaBar);
    phi = task.current_state;
    dx = task.grid.x(2)-task.grid.x(1);
    w2d = dx*dx;
    ndim = numel(size(task.grid.mesh.x));
    ddens = 2*imag(conj(phi).*(task.applyham(phi)-task.current_mu*phi)/(1+1i*task.gamma)); % - phi.*conj(task.applyham(phi)));
%     kinen = conj(phi).*task.grid.lap(phi);
    if(ndim==3)
        dz = task.grid.z(2)-task.grid.z(1);
        slicez = phi(:,:,task.grid.nz/2);
        densz = sum(abs(phi).^2,3)*dz;
        ddensz = sum((ddens),3)*dz;
%         kinen = sum((kinen),3)*dz;
        mulocal = sum(conj(phi).*task.applyham(phi),3)*(task.grid.z(2)-task.grid.z(1));
    else
        slicez = phi;
        densz = abs(phi).^2;
        ddensz = (ddens);
        mulocal = conj(phi).*task.applyham(phi);
    end

%     LLL = real(task.grid.inner(phi,task.grid.lz(phi)))./task.current_n;

    MU = real(sum(sum(mulocal))*w2d/task.current_n);
    
    ang1 = angle(sum(sum(slicez.*(angles3>(-pi + curangle)).*(angles3<-curangle))));
    dphase = angle(sum(sum(slicez.*exp(-1i*ang1).*(angles3<(-pi + curangle) | angles3>-curangle))));
    
    N1 = sum(sum(densz.*(angles3>(-pi + curangle)).*(angles3<-curangle)))*w2d;
    N2 = sum(sum(densz.*(angles3<(-pi + curangle) | angles3>-curangle)))*w2d;
    
%     K1 = sum(sum(kinen.*(angles3>(-pi + curangle)).*(angles3<-curangle)))*w2d;
%     K2 = sum(sum(kinen.*(angles3<(-pi + curangle) | angles3>-curangle)))*w2d;

    dN1 = sum(sum(ddensz.*(angles3>(-pi + curangle)).*(angles3<-curangle)))*w2d;
    dN2 = sum(sum(ddensz.*(angles3<(-pi + curangle) | angles3>-curangle)))*w2d;
    
    mu1 = real(sum(sum(mulocal.*(angles3>(-pi + curangle)).*(angles3<-curangle)))*w2d/N1);
    mu2 = real(sum(sum(mulocal.*(angles3<(-pi + curangle) | angles3>-curangle)))*w2d/N2);
    dmu = (mu1-mu2); 
    
%     maxn1 = max(max(densz.*(angles3>(-pi + curangle)).*(angles3<-curangle)));
%     maxn2 = max(max(densz.*(angles3<(-pi + curangle) | angles3>-curangle)));
    
    task.history.N1(task.current_iter) = N1;
    task.history.N2(task.current_iter) = N2;
    task.history.dN1(task.current_iter) = dN1;
    task.history.dN2(task.current_iter) = dN2;
%     task.history.K1(task.current_iter) = K1;
%     task.history.K2(task.current_iter) = K2;
%     task.history.maxn1(task.current_iter) = maxn1;
%     task.history.maxn2(task.current_iter) = maxn2;
    task.history.dmu(task.current_iter) = dmu;
    task.history.dphase(task.current_iter) = dphase;
    task.history.mu_int(task.current_iter) = MU;
%     save(sprintf('snapshots/slice_%05d',task.current_iter),'ddensz','slicez','mulocal','densz','curangle','N1','N2','dN1','dN2','dmu','dphase','MU');
%      imagesc(real(kinen)); 
%     caxis([-0.7 0.7]);
%     drawnow;
%     task.grid.imagesc(ddens);
%     caxis([-0.4 0.4]);
%     print('-dpng',sprintf('snapshots/ddens_%05d.png',task.current_iter), '-r150');
    txt = sprintf('z = %0.3f', (N1-N2)/task.Ntotal);
end
