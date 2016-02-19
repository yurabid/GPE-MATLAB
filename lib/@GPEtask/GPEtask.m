classdef GPEtask < handle
    %GPEtask - Solution of the Gross-Pitaevsii equation
    %   Detailed explanation goes here
    
    properties
        grid                   % grid object
        g                      % coupling coefficient
        omega = 0.0            % rotation speed
        gamma = 0.0            % dissipation constant
        decay_rate = 0         % decay rate
        Ntotal                 % initial total number of particles
        Vtrap                  % matrix of the trap potential
        Vtd                    % function handle to the time-dependent potential
        init_state             % initial state
        current_state          % current state in dynamics
        current_time = 0       % current time in dynamics
        current_iter = 0       % current iteration number in dynamics
        current_mu = 0         % current chemical potential in dynamics
        current_n              % current number of particles
        history                % history of current values
        user_callback          % user-defined callback function to process data after each step
    end
    
    methods
        function obj = GPEtask(grid,trappot)
            obj.grid = grid;
            obj.init_state = zeros(size(grid.mesh.x),'like',grid.mesh.x);
            if(isa(trappot,'double'))
                obj.Vtrap = trappot;
            elseif(isa(trappot,'function_handle'))
                obj.Vtrap = trappot(grid.mesh.x,grid.mesh.y,grid.mesh.z);
            end
            obj.dispstat('','init');
            obj.history = struct('mu',zeros(1,0,'gpuArray'),'n',zeros(1,0,'gpuArray'));
        end
        
        function v = getVtotal(obj,time)
            if(isa(obj.Vtd,'function_handle'))
                v = bsxfun(@plus,obj.Vtrap,obj.Vtd(obj.grid.mesh.x2,obj.grid.mesh.y2,time));
            else
                v = obj.Vtrap;
            end
        end
        
        function res = applyham(obj,phi)
            res = obj.grid.ifft(obj.grid.kk.*obj.grid.fft(phi)) + obj.getVtotal(obj.current_time).*phi + obj.g*abs(phi).^2.*phi;
        end
        
        function ext_callback(obj,phi,step,time,mu)
            if(exist('snapshots','file') ~= 7)
                mkdir('snapshots');
            end
            obj.current_state = phi;
            obj.current_time = time;
            obj.current_iter = step;
            obj.current_mu = mu;
            obj.history.mu(step) = mu;
            obj.current_n = obj.grid.norm(phi)^2;
            obj.history.n(step) = obj.current_n;

            if(isa(obj.user_callback,'function_handle'))
                obj.user_callback(obj);
            else
                ndim = numel(size(obj.grid.mesh.x));
                if(ndim==3)
                    slice = phi(:,:,obj.grid.nz/2);
                    densz = sum(abs(phi).^2,3)*(obj.grid.z(2)-obj.grid.z(1));
                else
                    slice = phi;
                    densz = abs(phi).^2;
                end
                save(sprintf('snapshots/slice_%05d',step),'slice','densz','time','mu');
            end
            ttime = toc;
            obj.dispstat(sprintf('Splitstep: iter - %u, mu - %0.3f, elapsed time - %0.3f seconds',step,mu,ttime));
        end
    end
    
end
