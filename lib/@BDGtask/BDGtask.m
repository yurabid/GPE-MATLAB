classdef BDGtask < handle
    % !!! Not properly tested, use with caution !!!
    % BDGtask - Solution of the Bogoliubov-de Gennes system of equations
    % Initialization:
    %     bdg = BDGtask(gpe_task), where gpe_task is a GPEtask object with
    %                              the ground-state solution initialized
    properties
       grid               % grid object
       ncomp = 1          % number of components
       g = 0
       kinen = {}
       nlin = []
       spectrum = []
       modes = []
       h0
       mtot
       mu = 0
       totsize = 0
       gridsize
       V
       phi
       ndim = 1
    end
    
    methods
        function obj = BDGtask(gpe_task,ncomp)
            if(nargin<2)
                ncomp = 1;
            end
            obj.grid = gpe_task.grid;
            obj.g = gpe_task.g;
            obj.ncomp = ncomp;
            obj.gridsize = size(obj.grid.mesh.x);
            obj.totsize = numel(obj.grid.mesh.x);
            obj.V = gpe_task.getVtotal(0);
            obj.phi = (gpe_task.init_state);
            obj.set_lap_cart(obj.grid.x,1);
            if(isprop(obj.grid,'y'))
                obj.set_lap_cart(obj.grid.y,2);
                obj.ndim=2;
            end
            if(isprop(obj.grid,'z'))
                obj.set_lap_cart(obj.grid.z,3);
                obj.ndim=3;
            end
            obj.set_nlin();
            obj.create_h0();
            obj.make_full_matr();            
        end
        
        function res=reshape_nto1(~,f)
            res=f(:);
        end
        function res=reshape_1ton(obj,f)
            res=reshape(f,obj.gridsize);
        end        
        function set_lap_cart(obj,grid,dim)
            obj.kinen{dim} = findiff1_arb(grid,2,0);
        end
        function set_lap_pol(obj,grid,dim,m)
            n = length(grid);
            obj.kinen{dim} = findiff1_arb(grid,2,0)+spdiags(1./grid(:),0,n,n)*...
                findiff1_arb(grid,1,0) - m^2*spdiags(1./grid(:).^2,0,n,n);
        end
        function set_nlin(obj)
            n = obj.totsize;
            obj.nlin = obj.g*spdiags(abs(obj.reshape_nto1(obj.phi)).^2,0,n,n);
        end
        function create_h0(obj)
            k = obj.kinen{1};
            if(obj.ndim>1)
                for i=2:obj.ndim
                    n1 = size(k,1);
                    n2 = size(obj.kinen{i},1);
%                     k = kron(speye(n2),k) + kron(obj.kinen{i},speye(n1));
                    k = kron(k,speye(n2)) + kron(speye(n1),obj.kinen{i});
                end
            end
            n = obj.totsize;
            obj.h0 = -0.5*k + spdiags(obj.reshape_nto1(obj.V),0,n,n);
        end
        function make_full_matr(obj)
            obj.get_mu2();
            n = obj.totsize;
            obj.mtot = [obj.h0 - obj.mu*speye(n) + 2*obj.nlin, obj.nlin;...
                -conj(obj.nlin), -obj.h0 - 2*obj.nlin + obj.mu*speye(n)];
        end 
        function mu=get_mu(obj)
            ph = obj.reshape_nto1(obj.phi);
            mu = real(ph'*(obj.h0 + obj.nlin)*ph./sum(abs(ph).^2));
            obj.mu = mu;
        end
        function mu=get_mu2(obj)
            mu = (eigs(obj.h0 + obj.nlin,1,0));
            obj.mu = mu;
        end        
        function [vv,sp]=solve(obj,nlev)
             [vv,dd] = eigs(obj.mtot,nlev,0);
             obj.spectrum = diag(dd);
             obj.modes = vv;
             sp = obj.spectrum;
        end
        function set_phi(obj,phi)
            obj.phi = phi;
            obj.set_nlin();
            obj.make_full_matr();
        end
        function set_v(obj,v)
            obj.V = v;
            obj.create_h0();
            obj.make_full_matr();
        end        
    end
    
end