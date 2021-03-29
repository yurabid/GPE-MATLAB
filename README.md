# GPE-MATLAB

Simple and extensible MATLAB solvers and tools for solving the Gross-Pitaevskii equation and related physical problems

### Installation: 
Add `lib` folder to MATLAB path (with subfolders).

### Some usage examples:
```matlab
% Initialization

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

% Stationary state calculation

tstep = 0.02; % initial time step for imaginary time evolution
acc = 1e-6; % desired accuracy
phi = task.groundstate_itp(tstep,acc); % phi will contain calculated stationary state

% Vortex dynamics calculation

xi = 1/sqrt(2*task.g*max(abs(phi(:))).^2); % healing length
x0 = 5; %vortex position
rr = sqrt((grid.mesh.x-x0).^2 + grid.mesh.y.^2);
angs = atan2(grid.mesh.x-x0,grid.mesh.y);
task.init_state = phi.*tanh(rr./xi).^1.*exp(1i*angs); % imprint a vortex on the stationary state
tstep = 0.004; % time step for calculation
steps_int = 50; % number of (internal) time steps in one (external) processing step
steps_ext = 500; % number of (external) processing steps
task.show_image = 1; % show the solution on every precessing step
task.solve_split(tstep,steps_int,steps_ext); % run the calculation

```

---------------------

Copyright (c) 2015-2021 Yuriy Bidasyuk (MIT License)

