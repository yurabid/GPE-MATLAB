function generate_particles_plummer(task, N, a, r0, rmax)
% Generates positions and velocities for N particles in a Plummer sphere.
% INPUT:
%   N: number of particles
%   a: Plummer scale length (default = 1)
%   r0: Coordinates of the center (default = [0,0,0])
%   rmax: Cut-off radius for the generated distribution

    task.bar_np = N;
    
    if nargin < 4
        r0 = [0,0,0]; % default center coords
    end
    if nargin < 3
        a = 1; % default Plummer scale length
    end
    % if nargin < 5
    rmin = 0; 
    % end
    if nargin < 5
        rmax = 100;
    end

    Mcoef = sqrt((task.bar_mass*(rmax^2+a^2)^(3/2)/rmax^3 + task.Ntotal)/2/pi);

    umin = (rmin/a)^3/(1+(rmin/a)^2)^(3/2);
    umax = (rmax/a)^3/(1+(rmax/a)^2)^(3/2);
    u = rand(N,1)*(umax-umin)+umin;
    r = a ./ sqrt(u.^(-2/3) - 1);
    
    theta = acos(2*rand(N,1) - 1);
    phi = 2*pi*rand(N,1);
    
    x = r .* sin(theta) .* cos(phi) + r0(1);
    y = r .* sin(theta) .* sin(phi) + r0(2);
    z = r .* cos(theta) + r0(3);
    task.bar_coords = [x, y, z];
    
    vel = zeros(N,3);
    for i = 1:N
        while true
            q = rand;
            p = rand;
            if p <= (1 - q^2)^(7/2) * q^2 / 0.1 % 0.1 is a safe upper bound for normalization
                break;
            end
        end
        v = q * Mcoef * (r(i)^2 + a^2)^(-1/4);
        theta_v = acos(2*rand - 1);
        phi_v = 2*pi*rand;
        vel(i,1) = v * sin(theta_v) * cos(phi_v);
        vel(i,2) = v * sin(theta_v) * sin(phi_v);
        vel(i,3) = v * cos(theta_v);
    end
    task.bar_vels = vel;
end
