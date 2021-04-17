%% Initialization 

units % load some common physical constants  
omegar = 150*2*pi; % trap frequencies (150Hz and 600Hz)
omegaz = 600*2*pi;
tcoef = 1/(omegar); % time scale
rscale = sqrt(hbar/mRB/omegar); % spatial scale

xmax = 15e-6/rscale; %  (15 micrometers)  grid boundaries [-xmax,xmax]
N=256; % number of grid points for coordinate grid
grid = grid2d(xmax,N,xmax,N); % create a two-dimensional NxN grid object

V = @(X,Y,Z) 0.5*(X.^2 + Y.^2); % function representing the external potential (must have 3 arguments)
task = GPEtask(grid,V); % initialize the GPE solver
task.Ntotal = 1.0e5; % wave function normalization (total number of atoms)
task.g = 4*pi*aRB/rscale/sqrt(2*pi/omegaz*omegar); % nonlinear interaction constant

%% Stationary state calculation

tstep = 0.02; % initial time step for imaginary time evolution
acc = 1e-6; % desired accuracy
phi = task.groundstate_itp(tstep,acc); % phi will contain calculated stationary state

%% Vortex dynamics calculation and detection

xi = 1/sqrt(2*task.g*max(abs(phi(:))).^2); % healing length
x0 = 5; %vortex position
rr = sqrt((grid.mesh.x-x0).^2 + grid.mesh.y.^2);
angs = atan2(grid.mesh.x-x0,grid.mesh.y);
task.init_state = task.init_state.*tanh(rr./xi).^1.*exp(1i*angs); % imprint a vortex at (5,0)
rr = sqrt((grid.mesh.x).^2 + (grid.mesh.y-x0).^2);
angs = atan2(grid.mesh.x,grid.mesh.y-x0);
task.init_state = task.init_state.*tanh(rr./xi).^1.*exp(1i*angs); % imprint another vortex at (0,5)
tstep = 0.004; % time step for calculation
steps_int = 50; % number of (internal) time steps in one (external) processing step
steps_ext = 500; % number of (external) processing steps
task.gamma=0.01; % non-zero damping constant turns solver to damped version of GPE
cms=[];
% simple way to implement post-processing is to run one external step at a
% time. This requires repetitive calls to solve_split increasing step number
for i=1:steps_ext
    task.solve_split(tstep,steps_int,i); % run the calculation
    imagesc(grid.x,grid.y,abs(task.current_state));hold on
    [cp,cm] = detect_core(task.current_state,grid.mesh.x,grid.mesh.y,1);
    plot(cm(:,1),cm(:,2),'+')
    hold off
    cms=[cms;cm]; % simplified way to collect core positions over time
    drawnow;
end
figure
plot(cms(:,1),cms(:,2),'.'); % plot trajectories of detected vortices 

%%
% More proper way to implement post-processing is to use a callback 
% functionality built into the solvers (also shows better performance)
task.current_iter = 0; % reset the iterations count to run the same dynamics again
task.user_callback = @post_process; % a function to be called on every processing step
task.history.coresx=[];
task.history.coresy=[];
task.solve_split(tstep,steps_int,steps_ext);
figure
plot(task.history.coresx,task.history.coresy,'.'); % plot trajectories of detected vortices 

function res = post_process(task)
    phi = task.current_state;
    [cp, cm] = detect_core(phi,task.grid.mesh.x,task.grid.mesh.y,1);
    task.grid.imagesc(abs(phi)); hold on;
    task.history.coresx = [task.history.coresx;cm(:,1)];
    task.history.coresy = [task.history.coresy;cm(:,2)];
    plot(cm(:,1),cm(:,2),'+');
    hold off;
    drawnow;
    res = '';
end