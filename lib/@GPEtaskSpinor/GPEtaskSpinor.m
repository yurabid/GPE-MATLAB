classdef GPEtaskSpinor < GPEtask
    %SGPEtask - Solution of the stochastic Gross-Pitaevsii equation
    
    properties
		coupling=0        % Inter-component coupling (Rabi frequency)
        ncomp=0           % Number of components      
    end
    
    methods
        function obj = GPEtaskSpinor(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
            if(isa(trappot,'cell'))
                obj.ncomp = numel(trappot);
                for i=1:obj.ncomp
                    if(isa(trappot{i},'function_handle'))
                        obj.Vtrap{i} = trappot{i}(grid.mesh.x,grid.mesh.y,grid.mesh.z);
                    end
                end
            end            
        end
        
        function v = getVtotal(obj,time)
            if(obj.ncomp == 0)
                obj.ncomp = size(obj.g,1);
            end
            if(time == 0 && numel(obj.V0)>1)
                    v = obj.V0;
            elseif(isa(obj.Vtd,'function_handle'))
                v = cell(1,obj.ncomp);
                for i=1:obj.ncomp
                    v{i} = bsxfun(@plus,obj.Vtrap{i},obj.Vtd(obj.grid.mesh.x2,obj.grid.mesh.y2,time));
                end
            elseif(isa(obj.Vtd,'cell'))
                v = cell(1,obj.ncomp);
                for i=1:obj.ncomp
                    v{i} = bsxfun(@plus,obj.Vtrap{i},obj.Vtd{i}(obj.grid.mesh.x2,obj.grid.mesh.y2,time));
                end                
            else
                v = obj.Vtrap;
            end
        end        
        function res = applyh0(obj,phi,j,time)
            if(nargin==3)
                time = obj.current_time;
            end
            v=obj.getVtotal(time);
            res = obj.grid.lap(phi) + v{j}.*phi;
            if(obj.omega ~= 0)
                res = res - obj.omega*obj.grid.lz(phi);
            end
        end        
        function res = applyham1c(obj,phi,j,time)
            if(nargin==3)
                time = obj.current_time;
            end            
            res = obj.getVtotal(time);
            res = res{j};
            for k=1:obj.ncomp
                res = res + obj.g(j,k)*abs(phi{k}).^2;
            end
            res = res.*phi{j} + obj.grid.lap(phi{j});
            for k=1:obj.ncomp
                if(j~=k)
                    res = res + obj.coupling*phi{k};
                end
            end
            if(obj.omega ~= 0)
                res = res - obj.omega*obj.grid.lz(phi{j});
            end            
        end
        
        function res = applyham(obj,phi,time)
            if(nargin==2)
                time = obj.current_time;
            end
            res = phi;
            for j =1:obj.ncomp
                res{j} = obj.applyham1c(phi,j,time);
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

