function generate_particles_hernquist(task, N, a)
% Generates positions and velocities for N particles in a Hernquist sphere.
% INPUT:
%   N: number of particles
%   a: Hernquist scale length
% OUTPUT:
%   pos: Nx3 matrix of positions (x,y,z)
%   vel: Nx3 matrix of velocities (vx,vy,vz)

    task.bar_np = N;
    M = (task.bar_mass)/4/pi;
    if nargin < 3
        a = 1; % default Hernquist scale length
    end

    % Generate positions 
    u = rand(N,1);
    r = a * sqrt(u) ./ (1 - sqrt(u));
    theta = acos(2*rand(N,1) - 1);
    phi = 2*pi*rand(N,1);
    x = r .* sin(theta) .* cos(phi);
    y = r .* sin(theta) .* sin(phi);
    z = r .* cos(theta);
    task.bar_coords  = [x, y, z];
    
    % Approximate velocity dispersion for Hernquist model (isotropic)
    sigma2 = M ./ (12*a) * (12*r.*(r + a).^3 ./ a^4) ./ ((r + a)/a).^4;
    sigma = sqrt(sigma2);
    
    % Sample 3D speed (Maxwellian approximation)
    v = sigma .* sqrt(-2 * log(rand(N,1)) - 2 * log(rand(N,1)) - 2 * log(rand(N,1)));
    
    % Sample direction (uniform on a sphere)
    theta_v = acos(2*rand(N,1) - 1);
    phi_v = 2*pi*rand(N,1);
    vel_x = v .* sin(theta_v) .* cos(phi_v);
    vel_y = v .* sin(theta_v) .* sin(phi_v);
    vel_z = v .* cos(theta_v);
    task.bar_vels = [vel_x, vel_y, vel_z];
end
