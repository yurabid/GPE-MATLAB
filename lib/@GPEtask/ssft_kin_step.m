function phi = ssft_kin_step(task,phi,dt)
    grid = task.grid;
    if(~isscalar(task.omega))
        phi = grid.ifft(task.kinop.*grid.fft(phi));
        lphi = phi;
        for ii = 1:task.n_crank
            lphi = phi + dt*task.omega.*grid.lz(lphi);
            lphi = 0.5*(phi+lphi);
        end
        phi = phi + dt*task.omega.*grid.lz(lphi);
        phi = grid.ifft(task.kinop.*grid.fft(phi));
    elseif(task.omega ~= 0)
        phi = grid.ifftx(task.kinop.*grid.fftx(phi));
        phi = grid.iffty(grid.ifftz(task.kinop2.*grid.fftz(grid.ffty(phi))));
        phi = grid.ifftx(task.kinop.*grid.fftx(phi));
    else
        phi = grid.ifft(task.kinop.*grid.fft(phi));
    end
end