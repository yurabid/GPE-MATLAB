function txt=post_process_p(task)
    global tshift rscale;
    phi = task.current_state;
    grid=task.grid;
    xm = 0.7e-6/rscale;
    curx = xm-bar_pos(task.current_time,xm,tshift);    
%     dx = task.grid.x(2)-task.grid.x(1);
%     w2d = dx*dx;
    ndim = numel(size(task.grid.mesh.x));
    phir = grid.sp2grid(phi);
%     imagesc(obj.grid.y,obj.grid.x,squeeze(abs(phir(:,:,end/2))));drawnow;
    if(ndim==3)
        slice = phir(:,:,end/2);
%         densz = sum(abs(phi).^2,3)*(task.grid.z(2)-task.grid.z(1));
%         mulocal = sum(conj(phi).*task.applyham(phi),3)*(task.grid.z(2)-task.grid.z(1));
    else
        slice = phir;
%         densz = abs(phi).^2;
%         mulocal = conj(phi).*task.applyham(phi);
    end
    
    N1 = grid.integrate(abs(grid.grid2sp(phir.*(task.grid.mesh.x<curx))).^2);
    N2 = grid.integrate(abs(grid.grid2sp(phir.*(task.grid.mesh.x>curx))).^2);
    NN1 = abs(phir).^2.*grid.wtot.*(task.grid.mesh.x<curx);
    NN2 = abs(phir).^2.*grid.wtot.*(task.grid.mesh.x>curx);
    ang1 = angle(sum(sum(slice.*(grid.mesh.x2<curx))));
    dphase = angle(sum(sum(slice.*exp(-1i*ang1).*(grid.mesh.x2>curx)))); 
    task.history.N1(task.current_iter) = N1;
    task.history.N2(task.current_iter) = N2;
    task.history.NN1(task.current_iter) = sum(NN1(:));
    task.history.NN2(task.current_iter) = sum(NN2(:));
    z1 = (N1-N2)/(N1+N2);
    z2 = (task.history.NN1(task.current_iter)-task.history.NN2(task.current_iter))/(task.history.NN1(task.current_iter)+task.history.NN2(task.current_iter));
    task.history.dphase(task.current_iter) = dphase;
%     LLL = -imag(task.grid.integrate(conj(phi).*(task.grid.mesh.x.*grid.derivy(phi) -...
%         task.grid.mesh.y.*task.grid.derivx(phi))))./task.current_n;

%     MU = real(sum(sum(mulocal))*w2d/task.current_n);
%     task.history.mu_int(task.current_iter) = MU;
%     imagesc(grid.y,grid.x,squeeze(abs(phir(:,:,end/2))));
%     h=pcolor(grid.z,grid.x,squeeze(abs(phir(:,end/2,:))));
%     set(h,'EdgeColor','none');
wpr = abs(phir).^2.*grid.mesh.wx.*grid.mesh.wy;
wpr1 = squeeze(sum(sum(wpr.*(grid.mesh.x<curx))));
wpr2 = squeeze(sum(sum(wpr.*(grid.mesh.x>curx))));
ncut=floor(length(grid.z)/4);
x = [ones(length(grid.z)-2*ncut,1) grid.z(ncut+1:end-ncut)'.^2];
b1=x\wpr1(ncut+1:end-ncut);
b2=x\wpr2(ncut+1:end-ncut);
task.history.b11(task.current_iter) = b1(1);
task.history.b12(task.current_iter) = b1(2);
task.history.b21(task.current_iter) = b2(1);
task.history.b22(task.current_iter) = b2(2);
% plot(grid.z,wpr1);
% hold on;
% plot(grid.z,wpr2);
% plot(grid.z(46:end-45),b1(1)+grid.z(46:end-45).^2*b1(2))
% plot(grid.z(46:end-45),b2(1)+grid.z(46:end-45).^2*b2(2))
%     [coresp, coresm] = detect_core(slice,grid.mesh.y2,grid.mesh.x2);
%     [coresp, coresm] = detect_core(phi.*((grid.mesh.x.^2+grid.mesh.y.^2)<16^2),task.grid.mesh.x,task.grid.mesh.y);
%     if(numel(coresp)>=4)
%          txt = 'NOTERM';
         
%     else
%         rads = coresp(:,1).^2+coresp(:,2).^2;
%         [v,i] = sort(rads);
%     %     plot(coresp(:,1),coresp(:,2),'MarkerSize',4,'Marker','x','LineStyle','none','Color',[1 1 1], 'LineWidth', 1);
%         plot(coresp(i(1:2),1),coresp(i(1:2),2),'MarkerSize',4,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1);
%         task.history.core1x(task.current_iter) = coresp(i(1),1);
%         task.history.core1y(task.current_iter) = coresp(i(1),2);
%         task.history.core2x(task.current_iter) = coresm(i(2),1);
%         task.history.core2y(task.current_iter) = coresm(i(2),2);        
    %     save(sprintf('snapshots/slice_%05d',task.current_iter),'slice','mulocal','densz','curx','N1','N2','dmu','dphase','MU');
%         task.grid.imagesc(abs(phi)); hold on;
if(mod(task.current_iter,20)==10)
    h=pcolor(grid.z,grid.y,squeeze(abs(phir(:,end/2,:))));
    set(h,'EdgeColor','none');
    drawnow;
end
if(mod(task.current_iter,20)==0)
    G=0;
%     phi2avg = 0;
%     phi4avg = 0;
    for i=1:length(task.snapshots)
        G = G+task.snapshots{i}*task.snapshots{i}';
%         ssc = phi*0;
%         ssc(grid.mask) = task.snapshots{i};
%         sscr = grid.sp2grid(ssc);
%         phi2avg = phi2avg + abs(sscr).^2;
%         phi4avg = phi4avg + abs(sscr).^4;
    end
    [task.current_c, task.current_nc] = eigs(G./length(task.snapshots),2);
%     phi2=phi*0;
%     phi2(grid.mask) = task.current_c;
    task.history.nc1(task.current_iter) = task.current_nc(1,1);
    task.history.nc2(task.current_iter) = task.current_nc(2,2);
%     phi2r = grid.sp2grid(phi2);
%     phi2avg = phi2avg./length(task.snapshots);
%     phi2avg = grid.sp2grid(phi2);
%     phi4avg = phi4avg./length(task.snapshots);
%     phi4avg = grid.sp2grid(phi2);
%     h=pcolor(grid.y,grid.x,squeeze(abs(phi2r(:,:,end/2))));
%     set(h,'EdgeColor','none');
%     plot(grid.x,abs(phi2r(:,end/2,end/2)).^2*task.current_nc);
%     hold on;
%     plot(grid.x,abs(phir(:,end/2,end/2)).^2);
%     plot(grid.x,abs(phi2avg(:,end/2,end/2)));
%     plot(grid.x,sqrt(abs(2*phi2avg(:,end/2,end/2).^2-phi4avg(:,end/2,end/2))));
%     hold off;
%     drawnow;
end
%     plot(grid.x,squeeze(abs(phir(:,end/2,end/2))));
%     ylim(gca,[0 25]);
%     end
%         hold off;
%         drawnow;
        if(numel(task.current_nc)>1)
            txt = sprintf('Z1 = %0.3f, Z2 = %0.3f, Nc = %0.3f, %0.3f', z1,z2, task.current_nc(1,1), task.current_nc(2,2));
        else
            txt = sprintf('Z1 = %0.3f, Z2 = %0.3f, Nc = %0.3f', z1,z2, task.current_nc);
        end
%     end

end
