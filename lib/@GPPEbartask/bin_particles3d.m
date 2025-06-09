function density = bin_particles3d(task, pos)
    M = task.grid.nx;
    N = size(pos, 1); % New number of particles
    % Find the indices of the lower corner of the cell containing each particle
    indices_low = floor(pos);
    
    % Calculate the fractional position within the cell
    delta = pos - indices_low;
    
    % Setup grid indices for the 8 corners of each cell
    x0 = indices_low(:,1) + 1;
    y0 = indices_low(:,2) + 1;
    z0 = indices_low(:,3) + 1;
    x1 = x0 + 1;
    y1 = y0 + 1;
    z1 = z0 + 1;
    
    % Calculate weights for trilinear interpolation
    wx0 = 1 - delta(:,1);
    wx1 = delta(:,1);
    wy0 = 1 - delta(:,2);
    wy1 = delta(:,2);
    wz0 = 1 - delta(:,3);
    wz1 = delta(:,3);
    
    % Initialize arrays to hold all indices and weights
    indicesAll = zeros(8*N, 1);
    weightsAll = zeros(8*N, 1);
    
    % Fill arrays for all 8 corners of each cell
    idx = 1:N;
    
    % Corner 000
    indicesAll(idx) = sub2ind([M, M, M], x0, y0, z0);
    weightsAll(idx) = wx0 .* wy0 .* wz0;
    idx = idx + N;
    
    % Corner 001
    indicesAll(idx) = sub2ind([M, M, M], x0, y0, z1);
    weightsAll(idx) = wx0 .* wy0 .* wz1;
    idx = idx + N;
    
    % Corner 010
    indicesAll(idx) = sub2ind([M, M, M], x0, y1, z0);
    weightsAll(idx) = wx0 .* wy1 .* wz0;
    idx = idx + N;
    
    % Corner 011
    indicesAll(idx) = sub2ind([M, M, M], x0, y1, z1);
    weightsAll(idx) = wx0 .* wy1 .* wz1;
    idx = idx + N;
    
    % Corner 100
    indicesAll(idx) = sub2ind([M, M, M], x1, y0, z0);
    weightsAll(idx) = wx1 .* wy0 .* wz0;
    idx = idx + N;
    
    % Corner 101
    indicesAll(idx) = sub2ind([M, M, M], x1, y0, z1);
    weightsAll(idx) = wx1 .* wy0 .* wz1;
    idx = idx + N;
    
    % Corner 110
    indicesAll(idx) = sub2ind([M, M, M], x1, y1, z0);
    weightsAll(idx) = wx1 .* wy1 .* wz0;
    idx = idx + N;
    
    % Corner 111
    indicesAll(idx) = sub2ind([M, M, M], x1, y1, z1);
    weightsAll(idx) = wx1 .* wy1 .* wz1;
    
    density_linear = accumarray(indicesAll, weightsAll, [M^3, 1]);
    density = reshape(density_linear, M, M, M);
end