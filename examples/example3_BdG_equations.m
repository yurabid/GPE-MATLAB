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

%% Calculation of the Bogoliubov-de Gennes collective excitations spectrum

bdg = BDGtask(task); % solver for BdG equations. It uses GPEtask as an input
[vv,sp]=bdg.solve(50); % calculate 50 lowest Bogoliubov modes

%%

% Solution contains both functions u and v from the solution. We plot only one of them
imagesc(real(bdg.reshape_1ton(vv(1:end/2,45)))); % show the shape of the mode #45. 
