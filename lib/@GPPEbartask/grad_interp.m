function [gradx, grady, gradz] = grad_interp(task, F, query_points)
    grid = task.grid;
    % [dFdx, dFdy, dFdz] = gradient(F, grid.dx, grid.dy, grid.dz);
    dFdx = grid.derivx(F);
    dFdy = grid.derivy(F);
    dFdz = grid.derivz(F);
    qx = query_points(:,1);
    qy = query_points(:,2);
    qz = query_points(:,3);
    gradx = interpn(grid.mesh.y, grid.mesh.x, grid.mesh.z, dFdx, qy, qx, qz, 'linear', 0);
    grady = interpn(grid.mesh.y, grid.mesh.x, grid.mesh.z, dFdy, qy, qx, qz, 'linear', 0);
    gradz = interpn(grid.mesh.y, grid.mesh.x, grid.mesh.z, dFdz, qy, qx, qz, 'linear', 0);
end