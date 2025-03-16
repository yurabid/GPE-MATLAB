classdef GPPEtask < GPEtask
    % GPPEtask - Solution of the Gross-Pitaevsii-Poisson system of equations
    % works only with 3D grids with equal dimensions
    
    properties
        Fi                       % stored value of gravitational potential
        use_QM_for_BC=true       % use Quadrupole momentum when calculating boundary conditions
        jacobi_niter=2           % number of Jacobi smoothing steps
        itp_fi_update_step=10    % ITP iterations between potential updates
        itp_use_fft_Fi=false     % ITP calculate grav. potential using FFT
        current_rcm              % Current center-of-mass position
        current_qm               % Current value of QM
        calculate_qmder=false    % Whether to calculate the third derivative of QM during evolution
        current_qmder            % Current value of the third derivative of QM
    end

    properties (Access=protected)
        ext_kernel
        jacobi_kernel
        jacobi_central_coef
        bcind
        bc_gridX
        bc_gridY
        bc_gridZ
        grid_inds
    end
    
    methods
        function obj = GPPEtask(grid,trappot)
            obj = obj@GPEtask(grid,trappot);

            obj.ext_kernel = zeros(3,3,3,'like',grid.x);
            obj.ext_kernel(:,:,1) = [1 2 1;2 4 2;1 2 1]/8;
            obj.ext_kernel(:,:,2) = [2 4 2;4 8 4;2 4 2]/8;
            obj.ext_kernel(:,:,3) = [1 2 1;2 4 2;1 2 1]/8;  

            obj.jacobi_kernel = zeros(3,3,3,'like',grid.x);

            % 19-point stncil
            % obj.jacobi_kernel(:,:,1) = [0 1 0;1 2 1;0 1 0]/24;
            % obj.jacobi_kernel(:,:,2) = [1 2 1;2 0 2;1 2 1]/24;
            % obj.jacobi_kernel(:,:,3) = [0 1 0;1 2 1;0 1 0]/24;
            % obj.jacobi_central_coef = 1/4;
            
            % 27-point stncil
            obj.jacobi_kernel(:,:,1) = [2 3 2;3 6 3;2 3 2]/88;
            obj.jacobi_kernel(:,:,2) = [3 6 3;6 0 6;3 6 3]/88;
            obj.jacobi_kernel(:,:,3) = [2 3 2;3 6 3;2 3 2]/88;
            obj.jacobi_central_coef = 26/88;

            obj.bcind=zeros(size(grid.mesh.x)+1);
            obj.bcind(1,:,:)=1;obj.bcind(end,:,:)=1;
            obj.bcind(:,1,:)=1;obj.bcind(:,end,:)=1;
            obj.bcind(:,:,1)=1;obj.bcind(:,:,end)=1;
            obj.bcind=obj.bcind(:)>0; 
            xg=[obj.grid.x, obj.grid.x(end)+obj.grid.x(2)-obj.grid.x(1)];
            [X,Y,Z]=meshgrid(xg,xg,xg);
            obj.bc_gridX=X(obj.bcind);
            obj.bc_gridY=Y(obj.bcind);
            obj.bc_gridZ=Z(obj.bcind);

            obj.grid_inds = 1:grid.nx;
        end
        
        function v = getVtotal(obj,time)
            v = obj.getVtotal@GPEtask(time) + obj.Fi(obj.grid_inds,obj.grid_inds,obj.grid_inds);
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
            tmp2 = obj.Fi;
            obj.Fi= obj.Fi*0.5;
            res = obj.get_energy@GPEtask(phi,time);
            obj.Fi = tmp2;
        end  
       
        function res=generate_pot(obj,phi)
            % generate monopole + quadrupole approximation
            % for gravitational potential
            dens = abs(phi).^2;
            Mtot = obj.grid.integrate(dens);
            xg=[obj.grid.x, obj.grid.x(end)+obj.grid.x(2)-obj.grid.x(1)];
            [X,Y,Z]=meshgrid(xg,xg,xg);
            res = obj.generate_pot_inner(dens,Mtot,X,Y,Z);
        end

        function res = generate_pot_inner(obj,dens,Mtot,Xg,Yg,Zg)
            % generate monopole + quadrupole approximation
            % for gravitational potential (inner function)
            rcm = obj.cm_coords(dens,Mtot);
            obj.current_rcm=rcm;
            X=Xg-rcm(1);
            Y=Yg-rcm(2);
            Z=Zg-rcm(3);
            rr = sqrt(X.^2+Y.^2+Z.^2);
            res=(-Mtot/(4*pi)./rr);
            if (obj.use_QM_for_BC)
                Q = obj.quad_mom(dens,rcm);
                obj.current_qm=Q;
                res = res - Q(1,1)/(8*pi)./rr.^5.*X.^2;
                res = res - Q(2,2)/(8*pi)./rr.^5.*Y.^2;
                res = res - Q(3,3)/(8*pi)./rr.^5.*Z.^2;
                res = res - Q(1,2)/(4*pi)./rr.^5.*X.*Y;
                res = res - Q(1,3)/(4*pi)./rr.^5.*X.*Z;
                res = res - Q(2,3)/(4*pi)./rr.^5.*Y.*Z;
            end
        end
        function set_pot(obj,phi,refine) 
            % Set gravitational potential based on arbitrary WF
            obj.Fi=obj.generate_pot(phi);
            if (nargin==3 && refine)
                h=obj.grid.x(2)-obj.grid.x(1);
                sz = size(obj.grid.mesh.x)+1;
                for i=1:10
                    [obj.Fi,~]=obj.V_cycle(obj.Fi,abs(phi).^2,h,sz(1));
                end
            end
        end        
        function set_pot_bc(obj,dens)
            % Set boundary conditions for gravitational potential
            obj.Fi(obj.bcind) = obj.generate_pot_inner(dens,obj.Ntotal,obj.bc_gridX,obj.bc_gridY,obj.bc_gridZ);
        end    
        function rcm = cm_coords(obj,dens,Mtot) 
            % calculate center-of-mass coordinates from density
            % distribution
            rcm=zeros(1,3,'like',dens);
            rcm(1) = obj.grid.integrate(dens.*obj.grid.mesh.x)/Mtot;
            rcm(2) = obj.grid.integrate(dens.*obj.grid.mesh.y)/Mtot;
            rcm(3) = obj.grid.integrate(dens.*obj.grid.mesh.z)/Mtot;            
            rcm=real(rcm);
        end        
        function Q = quad_mom(obj,dens,rcm)
            % Calculate quadrupole momentum (3x3 matrix) from density
            % distribution
            Q=zeros(3,3,'like',dens);
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
        
        function [Fi,res]=jacobi(obj,Fi,f,h,N)
            Fi_new = Fi;
            rhs = h^2*f(2:N-1,2:N-1,2:N-1)*obj.jacobi_central_coef;
            for u=1:obj.jacobi_niter
               Fi_new(2:N-1,2:N-1,2:N-1) = convn(Fi, obj.jacobi_kernel, 'valid') - rhs;
               if(u<obj.jacobi_niter)
                    Fi=Fi_new;
               end
            end
            res = (Fi-Fi_new)/(h^2*obj.jacobi_central_coef);
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
end

