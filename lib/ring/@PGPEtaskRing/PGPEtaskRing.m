classdef PGPEtaskRing < GPEtask
    %PGPEtask - Solution of the projected Gross-Pitaevsii equation
    
    properties
        T=0               % temperature
        snapshots={}      % snapshots for density matrix
        ss_maxlength=100  % maximal number of snapshots
        current_nc = 0
        current_c = 0
    end
    
    methods
        function obj = PGPEtaskRing(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
        end
        function res = applyham(obj,phi,time,phir)
            if(nargin==2)
                time = obj.current_time;
            end
            if(nargin<=3)
                 phir = obj.grid.sp2grid(phi);
            end            
            res = obj.grid.applyh0(phi) + obj.grid.grid2sp((obj.getVtotal(time) + obj.g*abs(phir.^2)).*phir);
        end
        
        function res = applyh0(obj,phi,time,phir)
           if(nargin==2)
                time = obj.current_time;
            end
            if(nargin<=3)
                 phir = obj.grid.sp2grid(phi);
            end  
            res = obj.grid.applyh0(phi) + obj.grid.grid2sp(obj.getVtotal(time).*phir);
        end        
    end
  
  methods (Access = protected) 
      
     
        function res=ext_callback(obj,phi,step,time,mu,n)
            if(exist('snapshots','file') ~= 7)
                mkdir('snapshots');
            end
            obj.snapshots{length(obj.snapshots)+1} = phi(obj.grid.mask);
            if(length(obj.snapshots)>obj.ss_maxlength)
                obj.snapshots = obj.snapshots(2:end);
            end
            phir = obj.grid.sp2grid(phi);
            V = obj.getVtotal(time);
            h0 = real(sum(obj.grid.to1d(conj(phi).*(obj.grid.applyh0(phi) + obj.grid.grid2sp(V.*phir)))));
            h1 = obj.grid.integrate_grid(obj.g*abs(phir).^4);
            mu = (h0 + h1);
            ee = (h0 + h1*0.5);            

            obj.current_state = phi;
            obj.current_time = time;
            obj.current_iter = step;
            obj.current_mu = mu;
            obj.history.mu(step) = mu;
            obj.current_n = n;
            obj.history.n(step) = n;
            obj.history.ee(step) = ee;
%             obj.history.nc(step) = obj.current_nc;
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
            obj.dispstat(sprintf(['Split-step: iter - %u, n - %0.3f, calc. time - %0.3f sec.; ',res_text],step,n,ttime));
            res = res_text;
        end
    end
    
end

