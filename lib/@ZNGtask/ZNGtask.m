classdef ZNGtask < GPEtask
    %ZNGtask - Solution of the Zaremba-Nikuni-Griffin equations
    
    properties
        T=0               % temperature
        init_state_nt         % initial state (thermal cloud)
        current_state_nt      % current state in dynamics (thermal cloud)
        current_nt        % current Nt
        current_nc        % current Nc
        grid_nt
        vtrap_nt
        n_test=1e6            % number of test particles
    end
    
    methods
        function obj = ZNGtask(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
            obj.grid_nt = grid3d(grid.x(end)*2,grid.nx*2,grid.y(end)*2,grid.ny*2,grid.z(end)*2,grid.nz*2);
            obj.vtrap_nt = trappot(obj.grid_nt.mesh.x,obj.grid_nt.mesh.y,obj.grid_nt.mesh.z);
        end
        
        function res=extend(obj,phi)
            nx = obj.grid.nx/2;
            ny = obj.grid.ny/2;
            nz = obj.grid.nz/2;
            res=zeros(4*nx,4*ny,4*nz,'like',obj.grid.x);
            res(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz) = phi;
        end
        
        function res=shrink(obj,phi)
            nx = obj.grid.nx/2;
            ny = obj.grid.ny/2;
            nz = obj.grid.nz/2;
            res=phi(nx+1:3*nx,ny+1:3*ny,nz+1:3*nz);
        end
        
        function res = applyham(obj,phi,nt,time)
            if(nargin==3)
                time = obj.current_time;
            end
            res = obj.getVtotal(time).*phi + obj.g.*(abs(phi).^2 + 2*obj.shrink(nt)).*phi;
            if(obj.omega ~= 0)
                res = res + obj.grid.lap(phi,obj.omega);
            else
                res = res + obj.grid.lap(phi);
            end
        end
    end
  
  methods (Access = protected) 
      
     
        function ext_callback(obj,phi,step,time,mu,n)
            if(exist('snapshots','file') ~= 7)
                mkdir('snapshots');
            end
            obj.current_state = phi;
            obj.current_time = time;
            obj.current_iter = step;
            obj.current_mu = mu;
            obj.history.mu(step) = mu;
            obj.current_n = n;
            obj.history.n(step) = n;
            res_text='';

            if(isa(obj.user_callback,'function_handle'))
                if(nargout(obj.user_callback) ~= 0)
                    res_text=obj.user_callback(obj);
                else
                    res_text='';
                    obj.user_callback(obj);
                end
            end
            ttime = toc;
            obj.dispstat(sprintf(['Split-step: iter - %u, mu - %0.3f, calc. time - %0.3f sec.; ',res_text],step,mu,ttime));
        end
    end
    
end

