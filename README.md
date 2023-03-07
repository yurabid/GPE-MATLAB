# GPE-MATLAB

Simple and extensible MATLAB solvers and tools for solving the Gross-Pitaevskii equation and related physical problems

### Installation: 
Add `lib` folder to MATLAB path (with subfolders).

### Structure of the package
There are two main types of objects in the project. In general, for any calculation you need one form each type.

1. Grid objects (class names start with 'grid') are used to define spatial grid for discretization of your solutions.
2. Solver objects (class names end with 'task') are used to calculare and analyze GPE and other related equations. These solvers can (mostly) work with any of the grid object type

Solver objects mostly inlude two main types of methods: groundstate methods (start with 'groundstate_') and dynamics methods (start with 'solve_').
Difference between methods of the same type is mostly in the algorithms used to solve the equation.

The list of available solvers for easier navigation:

- GPEtask -- Solver for a basic Gross-Pitaevskii equation. It is possible to add phenomenological damping to simulate energy dissipation or set up time-dependent external potentials.

- PGPEtask -- Projected Gross-Pitaevskii equation. Uses oscillator basis representation and restricted to three-dimensional simulations only. It is possible to add damping and noise to run Stochastic Projected Gross-Pitaevskii equation for finite-temperature simulations.

- GPEtask2comp -- same as GPEtask but for two-component condensate mixtures.

- GPEtaskSpinor -- same as GPEtask but for two-component coherently coupled condensates. The coherent coupling can be space-dependent and complex, allowing to simulate spin-orbit-coupled systems.


### Basic usage example:
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

### CUDA support:
Most generic solvers `GPEtask` and `GPEtaskSpinor` can also run on GPU (with MATLAB parallel computing toolbox).
To do this you need to initialize your grid with GPU arrays. For example

```matlab
xmax = 10;
N = 256;
x = gpuArray.linspace(-xmax,xmax,N);
grid = grid2d(x,x);
```
Other solvers should support such initialization as well, but they do not benefit too much from GPU parallelization.

---------------------

Copyright (c) 2015-2021 Yuriy Bidasyuk, [yurabid@gmail.com](mailto:yurabid@gmail.com) (MIT License)
