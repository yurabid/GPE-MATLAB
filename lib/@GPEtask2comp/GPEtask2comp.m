classdef GPEtask2comp < GPEtask
    %Solution of the two-component Gross-Pitaevsii equation for inceherent
    %mixture of two condensates
    
    properties 
        N1 % number of atoms in component 1
        N2 % number of atoms in component 2
        mu1 % chemical potential of component 1
        mu2 % chemical potential of component 2
        M1=1 % atom mass of component 1
        M2=1 % atom mass of component 2
        ncomp=2
    end
    
    methods
        function obj = GPEtask2comp(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
            if(isa(trappot,'cell'))
                ncomp = numel(trappot);
                for i=1:ncomp
                    if(isa(trappot{i},'function_handle'))
                        obj.Vtrap{i} = trappot{i}(grid.mesh.x,grid.mesh.y,grid.mesh.z);
                    end
                end
            end            
        end
        
        function v = getVtotal(obj,time)
            if(isempty(obj.Vtd))
                v = obj.Vtrap;
            elseif(time == 0 && numel(obj.V0)>1)
                v = obj.V0;
            elseif(isa(obj.Vtd,'function_handle'))
                v = {obj.Vtrap{1} + obj.Vtd(obj,time), obj.Vtrap{2} + obj.Vtd(obj,time)};
            elseif(isa(obj.Vtd,'cell'))
                v = {obj.Vtrap{1} + obj.Vtd{1}(obj,time), obj.Vtrap{2} + obj.Vtd{2}(obj,time)};
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
        
        function res = inner(obj,left,right)
            res = 0;
            for j =1:obj.ncomp
                res = res + obj.grid.inner(left{j},right{j});
            end
        end
        
        function res = mul(~,phi,a)
            res = phi;
            res{1} = res{1}*a;
            res{2} = res{2}*a;
        end
        
        function res = get_energy(obj,phi,time)
            if(nargin<3)
                time = obj.current_time;
            end
            if(nargin<2)
                phi = obj.current_state;
            end            
            tmp = obj.g;
            obj.g = 0.5*obj.g;
            res = real(obj.inner(phi,obj.applyham(phi,time)));
            obj.g=tmp;
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
            if(obj.show_image>0)
                 hold off; obj.grid.imagesc(abs(phi{obj.show_image})); drawnow;
            end
            ttime = toc;
            obj.dispstat(sprintf(['Split-step: iter - %u, mu - %0.3f, calc. time - %0.3f sec.; ',res_text],step,mu,ttime));
        end
    end
    
end

