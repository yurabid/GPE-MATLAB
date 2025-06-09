classdef GPPEbartask < GPPEtask
    % GPPEbartask - Solution of the Gross-Pitaevsii-Poisson system of
    % equations + self-consistent particle-mesh calculation of the gas 
    % of classical particles gravitationally interacting with the condensate.  
    % works only with 3D grids with equal dimensions
    
    properties
        bar_coords = [0,0,0]
        bar_vels = [0,0,0]
        bar_dens = 0
        bar_Fi = 0
        bar_mass = 0
        bar_np = 1
        bar_amp = 0
        bar_rad = 0
        solve_bec_dynamic=true
        solve_bar_dynamic=true
    end
    
    methods
        function obj = GPPEbartask(grid,trappot)
            obj = obj@GPPEtask(grid,trappot);
            obj.bar_dens=zeros(size(grid.mesh.x),'like',grid.mesh.x);
        end

        function generate_particles_disc(obj, np, rmax, zmax)
            obj.bar_np = np;
            rs = (rand(np,1)).^(0.55)*rmax;
            zs = rand(np,1)*2*zmax-zmax;
            phis = rand(np,1)*2*pi;

            obj.bar_coords = [rs.*cos(phis), rs.*sin(phis), zs];
            obj.bar_vels = zeros(np,3);
            obj.bar_dens=obj.particle_cloud_density()/obj.grid.weight*obj.bar_mass/obj.bar_np;
            rgrid = (0:99)*rmax/100;
            vgrid = rgrid*0;
            rr=sqrt(obj.grid.mesh.x.^2 + obj.grid.mesh.y.^2);
            for i=1:100
                mcur = obj.grid.integrate(obj.current_state(rr<rgrid(i))+obj.bar_dens(rr<rgrid(i)));
                vgrid(i) = real(sqrt(mcur/4/pi/rgrid(i)));
            end
            obj.bar_vels = (interp1(rgrid,vgrid,rs)).*[-sin(phis)- rand(np,1)*0.4 + 0.2,...
                cos(phis)- rand(np,1)*0.4 + 0.2, rand(np,1)*0.4 - 0.2];
            % obj.bar_vels = (2+2*rand(np,1)).*[-sin(phis), cos(phis), zeros(np,1)];
            % obj.bar_vels = (sqrt(obj.Ntotal/4/pi + obj.bar_mass/8/pi) + rand(np,1)*0.4-0.2).*[-sin(phis), cos(phis), zeros(np,1)].*rs./(rs.^1.5+2);
        end

        function density = bin_particles1d(task, pos, weights)
            if nargin < 3
                weights = ones(size(pos));
            end
            indices_low = floor(pos);
            delta = pos - indices_low;
            x0 = indices_low(:) + 1;
            M = max(indices_low) + 2;
            density = accumarray([x0; x0 + 1], [(1 - delta(:)).*weights(:); delta(:).*weights(:)], [M, 1]);
        end

        function [rdist, rgrid] = bar_density_vs_radius(task)
            rstep = task.grid.dx;
            rvals = sqrt(task.bar_coords(:,1).^2+task.bar_coords(:,2).^2+task.bar_coords(:,3).^2);
            rdist = task.bin_particles1d(rvals./rstep)/task.bar_np*task.bar_mass/rstep;
            rgrid = (0:numel(rdist)-1)*rstep;
        end

        function [vel_avg, sigma, rgrid] = bar_velocity_vs_radius(task)
            rstep = task.grid.dx;
            pind = sqrt(task.bar_coords(:,1).^2+task.bar_coords(:,2).^2+task.bar_coords(:,3).^2)./rstep;
            indices_low = floor(pind);
            M = max(indices_low) + 1;
            % n_per_bin = groupcounts(indices_low);
            velocity = sqrt(task.bar_vels(:,1).^2+task.bar_vels(:,2).^2+task.bar_vels(:,3).^2);
            % vel_avg = accumarray(indices_low+1,velocity, [M, 1])./n_per_bin;
            rgrid = (0:M-1)*rstep;
            vel_avg = zeros(M,1);
            sigma = zeros(M,1);
            for i=1:M
                vel_in_bin = velocity(indices_low==i-1);
                n_per_bin = numel(vel_in_bin);
                if n_per_bin > 0
                    vel_avg(i) = sum(vel_in_bin)/n_per_bin;
                    sigma(i) = sqrt(sum((vel_in_bin-vel_avg(i)).^2)/n_per_bin)*sqrt(3);
                end
            end
        end

        function res = get_energy(obj,varargin)
            tmp2 = obj.Fi;
            obj.Fi= obj.Fi(obj.grid_inds,obj.grid_inds,obj.grid_inds)*0.5+obj.bar_Fi*0.5;
            if obj.bar_mass > 0
                m0 = obj.bar_mass/obj.bar_np;
                res = obj.get_energy@GPEtask(varargin{:});
                res = res + obj.grid.integrate(obj.Fi.*obj.bar_dens);
                res = res + sum(obj.bar_vels(:).^2)/2*m0;
            end
            obj.Fi = tmp2;
        end
        
        function bdens = bar_fun(obj,bc)
            bdens = obj.bar_amp*exp(-((obj.grid.mesh.x-bc(1)).^2+...
                (obj.grid.mesh.y-bc(2)).^2+(obj.grid.mesh.z-bc(3)).^2)/obj.bar_rad^2);
        end
        
        function res=generate_bar_dens(obj,bc)
            if(nargin<2)
                bc=obj.bar_coords;
            end
            np=size(bc,1);
            res=0;
            for i=1:np
                res = res + obj.bar_fun(bc(i,:));
            end
            obj.bar_dens=res;
        end

        function set_pot_bar_an(obj,rcm)
            X = obj.grid.mesh.x-rcm(1);
            Y = obj.grid.mesh.y-rcm(2);
            Z = obj.grid.mesh.z-rcm(3);
            R = sqrt(X.^2+Y.^2+Z.^2);
            obj.bar_Fi = -obj.bar_mass/4/pi*erf(R/obj.bar_rad)./R;
        end  

        function set_bar_mass_amp(obj)
            if obj.bar_amp ~= 0 && obj.bar_mass == 0
                obj.bar_mass = obj.bar_amp*pi^1.5*obj.bar_rad^3;
            elseif obj.bar_amp == 0 && obj.bar_mass ~= 0
                obj.bar_amp = obj.bar_mass/(pi^1.5*obj.bar_rad^3);
            end
        end      
    end
end
