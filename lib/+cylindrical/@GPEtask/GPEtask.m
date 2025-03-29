classdef GPEtask < GPEtask
    %GPEtask - Solution of the Gross-Pitaevsii equation
    
    properties

    end
    
    methods
        function obj = GPEtask(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
            if(isa(trappot,'function_handle'))
                obj.Vtrap = trappot(grid.mesh.r,grid.mesh.phi,grid.mesh.z);
            else
                obj.Vtrap = trappot;
            end
        end
        
        function set_kinop(task,dt)
            grid = task.grid;
            task.kinop2 = exp(-(grid.kr.^2+grid.kz.^2)/2*dt);
            task.kinop = exp(-((grid.kphi.^2-1/4)./grid.mesh.r.^2)/4*dt);
        end

        function phi = ssft_kin_step(task,phi,~)
            grid = task.grid;
            phi = phi.*sqrt(grid.mesh.r);
            phi = grid.ifftphi(task.kinop.*grid.fftphi(phi));
            phi = grid.ifftr(grid.ifftz(task.kinop2.*grid.fftz(grid.fftr(phi))));
            phi = grid.ifftphi(task.kinop.*grid.fftphi(phi));
            phi = phi./sqrt(grid.mesh.r);
        end
       
        function res = applyham(obj,phi,time)
            arguments
                obj GPEtask
                phi = obj.current_state
                time = obj.current_time
            end
            res = obj.getVtotal(time).*phi + obj.g*abs(phi).^2.*phi;
            res = res + obj.grid.lap(phi,obj.omega);
        end
        
        function res = applyh0(obj,phi,time)
            arguments
                obj GPEtask
                phi = obj.current_state
                time = obj.current_time
            end
            res = obj.getVtotal(time).*phi;
            res = res + obj.grid.lap(phi,obj.omega);
        end
       
        
    end
  
end

