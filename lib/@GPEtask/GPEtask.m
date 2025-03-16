classdef GPEtask < handle
    %GPEtask - Solution of the Gross-Pitaevsii equation

    properties
        grid               % grid object
        g                  % coupling coefficient
        omega = 0.0        % rotation speed
        n_crank = 3        % number of Crank-Nicolson iterations for L
        n_recalc = 10      % number of iterations to recalc potential and chem.pot.
        gamma = 0.0        % dissipation constant
        decay_rate = 0     % 1/e decay time (0 for no decay)
        Ntotal = 0         % initial total number of particles
        mu_init = 0        % initial chemical potential
        Vtrap              % matrix of the static trap potential
        Vtd                % function handle to the time-dependent potential
        V0 = 0             % INTERNAL: total potential at t=0
        Vcurrent = 0       % matrix of the current total trap potential
        init_state         % initial state
        current_state      % current state in dynamics
        current_time = 0   % current time in dynamics
        current_iter = 0   % current iteration number in dynamics
        current_mu = 0     % current chemical potential in dynamics
        current_n          % current number of particles
        history            % history of current values
        user_callback      % user-defined callback function to process data after each step
        itp_max_iter=5000  % ITP maximum iterations
        fftw_set_wisdom=true
        itp_adjust_stepsize=false
        % real time propagation properties
        dt                 % time step for history arrays
        totalTime          % total evolution time
        show_image = 0     % show density image on each time step
    end
    properties (Access=protected)
        kinop
        kinop2
        nlinop
    end
    methods
        function obj = GPEtask(grid,trappot)
            obj.grid = grid;
            obj.init_state = zeros(size(grid.mesh.x),'like',grid.x);
            if(isa(trappot,'function_handle'))
                obj.Vtrap = trappot(grid.mesh.x,grid.mesh.y,grid.mesh.z);
            else
                obj.Vtrap = trappot;
            end
            rng('shuffle');
            if(isa(grid.mesh.x,'gpuArray'))
                parallel.gpu.rng('shuffle');
            end
            obj.dispstat('','init');
            obj.history = struct('mu',zeros(1,0,'like',grid.x),'n',zeros(1,0,'like',grid.x));
            if (obj.fftw_set_wisdom)
                fftw('dwisdom',[]);
                fftw('planner','measure');
                fftn(grid.mesh.x);
            end
        end

        function v = getVtotal(obj,time)
            if(time == 0 && numel(obj.V0)>1)
                v = obj.V0;
            elseif(isa(obj.Vtd,'function_handle'))
%                 v = bsxfun(@plus,obj.Vtrap,obj.Vtd(obj.grid.mesh.x2,obj.grid.mesh.y2,time));
                v = obj.Vtrap + obj.Vtd(obj,time);
            else
                v = obj.Vtrap;
            end
        end

        function res = applyham(obj,phi,time)
            if(nargin==2)
                time = obj.current_time;
            end
            res = obj.applyh0(phi,time) + obj.g.*abs(phi).^2.*phi;
        end

        function res = applyh0(obj,phi,time)
            if(nargin==2)
                time = obj.current_time;
            end
            res = obj.getVtotal(time).*phi;
            if(isscalar(obj.omega) && obj.omega ~= 0)
                res = res + obj.grid.lap(phi,obj.omega);
            else
                res = res + obj.grid.lap(phi);
            end
            if(~isscalar(obj.omega))
                res = res - obj.omega.*obj.grid.lz(phi);
            end 
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
            res = real(obj.grid.inner(phi,obj.applyham(phi,time)));
            obj.g=tmp;
        end

  end

  methods (Access = protected)

      dispstat(obj,TXT,varargin);

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
            if(obj.show_image>0)
                hold off; obj.grid.imagesc(abs(phi)); drawnow;
            end
            ttime = toc;
            obj.dispstat(sprintf(['Split-step: iter - %u, mu - %0.3f, calc. time - %0.3f sec.; ',res_text],step,mu,ttime));
            res = res_text;
        end
    end

end

