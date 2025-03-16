function set_kinop(task,dt)
    grid = task.grid;
    if(~isscalar(task.omega))
        task.kinop = exp(-grid.kk*dt/2);
    elseif(task.omega ~= 0)
        task.kinop = exp(-(grid.kx.^2/2-grid.kx.*grid.mesh.y*task.omega)*dt/2);
        task.kinop2 = exp(-(grid.ky.^2/2+grid.kz.^2/2+grid.ky.*grid.mesh.x*task.omega)*dt);
    else
        task.kinop = exp(-grid.kk*dt);
    end  
end