classdef GPPEtask < GPEtask
    % GPPEtask - Solution of the Gross-Pitaevsii-Poisson system of equations
    % works only with 3D grids with equal dimensions
    
    properties
        alpha=0
        nl_amp=1
        Fi
        ext_kernel
        laps={}
        trans={}
        sizes={}
        bcind
        current_qm
        current_qmder
        current_rcm
        bar_dens
        bar_Fi
    end
    
    methods
        function obj = GPPEtask(grid,trappot)
            obj = obj@GPEtask(grid,trappot);
            obj.ext_kernel = zeros(3,3,'like',grid.x);
            obj.ext_kernel(:,:,1) = [1 2 1;2 4 2;1 2 1]/8;
            obj.ext_kernel(:,:,2) = [2 4 2;4 8 4;2 4 2]/8;
            obj.ext_kernel(:,:,3) = [1 2 1;2 4 2;1 2 1]/8;  
            obj.bcind=zeros(size(grid.mesh.x)+1);
            obj.bcind(1,:,:)=1;obj.bcind(end,:,:)=1;
            obj.bcind(:,1,:)=1;obj.bcind(:,end,:)=1;
            obj.bcind(:,:,1)=1;obj.bcind(:,:,end)=1;
            obj.bcind=obj.bcind(:)>0;  
            obj.bar_dens=zeros(size(grid.mesh.x),'like',grid.mesh.x);
        end
        
        function res = applyham(obj,phi,time)
            if(nargin==2)
                time = obj.current_time;
            end
            res = (obj.getVtotal(time)+obj.Fi(1:end-1,1:end-1,1:end-1)).*phi + obj.g.*abs(phi).^2.*phi;
            res = res + obj.grid.lap(phi);
        end

        function res = applyh0(obj,phi,time)
            if(nargin==2)
                time = obj.current_time;
            end
            res = (obj.getVtotal(time)+obj.Fi(1:end-1,1:end-1,1:end-1)).*phi;
            res = res + obj.grid.lap(phi);
        end
        function res = get_kin_energy(obj,phi)
            res = real(obj.grid.inner(phi,obj.grid.lap(phi)));
        end        
        function res = get_energy(obj,phi,time)
            if(nargin<3)
                time = obj.current_time;
            end
            if(nargin<2)
                phi = obj.current_state;
            end
            tmp = obj.g;
            tmp2 = obj.Fi;
            obj.g = 0.5*obj.g;
            obj.Fi= obj.Fi*0.5;
            res = real(obj.grid.inner(phi,obj.applyham(phi,time)));
            obj.g=tmp;
            obj.Fi = tmp2;
        end        
        function res=generate_pot(obj,phi)
%             phi = sqrt(obj.Ntotal)*obj.grid.normalize(phi);
            dens = abs(phi).^2+obj.bar_dens;
            Mtot = obj.grid.integrate(dens);
            rcm = obj.cm_coords(dens,Mtot);
            obj.current_rcm=rcm;
            xg=[obj.grid.x, obj.grid.x(end)+obj.grid.x(2)-obj.grid.x(1)];
            [X,Y,Z]=meshgrid(xg-rcm(1),xg-rcm(2),xg-rcm(3));
            rr = sqrt(X.^2+Y.^2+Z.^2);
            res=(-Mtot/(4*pi)./rr);
 
            Q = obj.quad_mom(dens,rcm);
            obj.current_qm=Q;
            res = res - Q(1,1)/(8*pi)./rr.^5.*X.^2;
            res = res - Q(2,2)/(8*pi)./rr.^5.*Y.^2;
            res = res - Q(3,3)/(8*pi)./rr.^5.*Z.^2;
            res = res - Q(1,2)/(4*pi)./rr.^5.*X.*Y;
            res = res - Q(1,3)/(4*pi)./rr.^5.*X.*Z;
            res = res - Q(2,3)/(4*pi)./rr.^5.*Y.*Z;
        end
        function set_pot(obj,phi,refine)            
            obj.Fi=obj.generate_pot(phi);
            if (nargin==3 && refine)
                h=obj.grid.x(2)-obj.grid.x(1);
                sz = size(obj.grid.mesh.x)+1;
                for i=1:5
                    [obj.Fi,~]=obj.V_cycle(obj.Fi,abs(phi).^2,h,sz(1));
                end
            end
        end        
        function set_pot_bc(obj,phi)
            tmp = obj.generate_pot(phi);
            obj.Fi(obj.bcind) = tmp(obj.bcind);
        end    
        function rcm = cm_coords(obj,dens,Mtot)
            rcm=zeros(1,3,'like',dens);
%             dens = abs(phi).^2+task.bar_dens;
            rcm(1) = obj.grid.integrate(dens.*obj.grid.mesh.x)/Mtot;
            rcm(2) = obj.grid.integrate(dens.*obj.grid.mesh.y)/Mtot;
            rcm(3) = obj.grid.integrate(dens.*obj.grid.mesh.z)/Mtot;            
            rcm=real(rcm);
        end        
        function Q = quad_mom(obj,dens,rcm)
            Q=zeros(3,3,'like',dens);
%             dens = abs(phi).^2+task.bar_dens;
            X = (obj.grid.mesh.x-rcm(1));
            Y = (obj.grid.mesh.y-rcm(2));
            Z = (obj.grid.mesh.z-rcm(3));
            rr = (X.^2+Y.^2+Z.^2);
            Q(1,2) = obj.grid.integrate(dens.*(3*X.*Y));
            Q(1,3) = obj.grid.integrate(dens.*(3*X.*Z));
            Q(2,3) = obj.grid.integrate(dens.*(3*Z.*Y));
            
            Q(2,1) = Q(1,2); Q(3,1) = Q(1,3); Q(3,2) = Q(2,3);
            
            Q(1,1) = obj.grid.integrate(dens.*(3*X.^2-rr));
            Q(2,2) = obj.grid.integrate(dens.*(3*Y.^2-rr));
            Q(3,3) = obj.grid.integrate(dens.*(3*Z.^2-rr));
            Q=real(Q);
        end
      
        function [phi,r] = V_cycle(obj,phi,f,h,N,Nmin)
            if nargin<6
                Nmin=3;
            end
            [phi,r] = obj.jacobi(phi,f,h,N);
            rhs = obj.restrict(r);
            sz = size(rhs);
            eps = zeros(sz,'like',rhs);
            sz=sz(1);
            if sz<=Nmin
                [eps,~] = obj.jacobi(eps,rhs,2*h,sz);
            else
                [eps,~] = obj.V_cycle(eps,rhs,2*h,sz,Nmin);
            end
            phi = phi + obj.extend(eps);
            [phi,r] = obj.jacobi(phi,f,h,N);
        end    
        
        function [Fi,res]=jacobi(~,Fi,f,h,N)
            Fi_new = Fi;
            niter=3;% number of smoothing steps, 2 seems to be enough
            for u=1:niter
               Fi_new(2:N-1,2:N-1,2:N-1)=(...
                   2*(  Fi(3:N,2:N-1,2:N-1) + Fi(1:N-2,2:N-1,2:N-1)...
                      + Fi(2:N-1,3:N,2:N-1) + Fi(2:N-1,1:N-2,2:N-1)...
                      + Fi(2:N-1,2:N-1,3:N) + Fi(2:N-1,2:N-1,1:N-2))...
                   - 6*h^2*f(2:N-1,2:N-1,2:N-1)...
                   + (Fi(3:N,3:N,2:N-1)    + Fi(1:N-2,3:N,2:N-1)...
                   + Fi(3:N,1:N-2,2:N-1)  + Fi(1:N-2,1:N-2,2:N-1)...
                   + Fi(2:N-1,3:N,3:N)    + Fi(2:N-1,1:N-2,3:N)...
                   + Fi(2:N-1,3:N,1:N-2)  + Fi(2:N-1,1:N-2,1:N-2)...
                   + Fi(3:N,2:N-1,3:N)    + Fi(3:N,2:N-1,1:N-2)...
                   + Fi(1:N-2,2:N-1,3:N)  + Fi(1:N-2,2:N-1,1:N-2))...                  
                   )/24; 
               if(u<niter)
                    Fi=Fi_new;
               end
            end
            res = (Fi-Fi_new)/(h^2/4);
        end    

        function [Fi,res]=jacobi2(~,Fi,f,h,N)
            Fi_new = Fi;
            niter=3;% number of smoothing steps, 2 seems to be enough
            for u=1:niter
               Fi_new(2:N-1,2:N-1,2:N-1)=(...
                   6*(  Fi(3:N,2:N-1,2:N-1) + Fi(1:N-2,2:N-1,2:N-1)...
                      + Fi(2:N-1,3:N,2:N-1) + Fi(2:N-1,1:N-2,2:N-1)...
                      + Fi(2:N-1,2:N-1,3:N) + Fi(2:N-1,2:N-1,1:N-2))...
                   - 26*h^2*f(2:N-1,2:N-1,2:N-1)...
                   + 3*(Fi(3:N,3:N,2:N-1)    + Fi(1:N-2,3:N,2:N-1)...
                   + Fi(3:N,1:N-2,2:N-1)  + Fi(1:N-2,1:N-2,2:N-1)...
                   + Fi(2:N-1,3:N,3:N)    + Fi(2:N-1,1:N-2,3:N)...
                   + Fi(2:N-1,3:N,1:N-2)  + Fi(2:N-1,1:N-2,1:N-2)...
                   + Fi(3:N,2:N-1,3:N)    + Fi(3:N,2:N-1,1:N-2)...
                   + Fi(1:N-2,2:N-1,3:N)  + Fi(1:N-2,2:N-1,1:N-2))...
                   + 2*(Fi(3:N,3:N,3:N)    + Fi(1:N-2,3:N,3:N)...
                   + Fi(3:N,3:N,1:N-2)  + Fi(1:N-2,3:N,1:N-2)...
                   + Fi(3:N,1:N-2,3:N)    + Fi(1:N-2,1:N-2,3:N)...
                   + Fi(3:N,1:N-2,1:N-2)  + Fi(1:N-2,1:N-2,1:N-2))...                   
                   )/88;
                
               if(u<niter)
                    Fi=Fi_new;
               end
            end
            res = (Fi-Fi_new)/(h^2/88*26);
        end          
        
        function phi=restrict(~,phi)
            phi=phi(1:2:end,1:2:end,1:2:end);
        end

        function res=extend(obj,phi)
            res = zeros(size(phi)*2-1,'like',phi); 
            res(1:2:end,1:2:end,1:2:end) = phi;
            res = convn(res,obj.ext_kernel,'same');
        end  
    end
  
  methods (Access = protected)         
        function res=ext_callback(obj,phi,step,time,mu,n)
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
            res = res_text;
        end
    end
    
end

