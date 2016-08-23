classdef PGPEtask < GPEtask
    %PGPEtask - Solution of the projected Gross-Pitaevsii equation
    
    properties
        T=0               % temperature
    end
    
    methods
        function obj = PGPEtask(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
        end
        function res = applyham(obj,phi,time,phir)
            if(nargin==2)
                time = obj.current_time;
            end
            if(nargin<=3)
                 phir = obj.grid.sp2grid(phi);
            end            
            res = obj.grid.etot.*phi + obj.grid.grid2sp((obj.getVtotal(time) + obj.g*abs(phir.^2)).*phir);
        end
        
        function res = applyh0(obj,phi,time,phir)
           if(nargin==2)
                time = obj.current_time;
            end
            if(nargin<=3)
                 phir = obj.grid.sp2grid(phi);
            end  
            res = obj.grid.etot.*phi + obj.grid.grid2sp(obj.getVtotal(time).*phir);
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
            phir = obj.grid.sp2grid(phi);
            imagesc(obj.grid.y,obj.grid.x,squeeze(abs(phir(:,:,end/2))));drawnow;
        end
    end
    
end

