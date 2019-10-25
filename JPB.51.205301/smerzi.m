% global tcoef Ub omegaBar initphase ub1 r0 r1 r01;
units
global tcoef Ub rscale vv;

% Smerzi PhysRevLett.84.4521
aRB = 58.19e-10;
omegascale = 50*2*pi;
omz = 17.68/50;
tcoef = 1/(omegascale); %time scale
rscale = sqrt(hbar/mRB/omegascale); % length scale
vv = 0;
Ub=650/omegascale*2*pi; % barrier height
totaltime = 1/tcoef;
NN = 5.0e4; % number of particles
g = 4*pi*aRB/rscale; % nonlinear coupling strength
Vfun = @(X,Y,Z) 0.5*(omz^2*Z.^2 + omz^2*Y.^2 + X.^2); % static part of the potential

% Initialize the grid
xmax = 12;
zmax = 12;
N=64;
Nz=64;
grid = grid3dgpu(xmax,N,xmax,N,zmax,Nz);
% Initialize the solver
task = GPEtask(grid,Vfun);
task.g = g;
task.Ntotal = NN;
task.Vtd = @TDPot; % time dependent part of the potential
task.user_callback = @post_process_sample;
% Calculate the ground state
[phi, mu1] = task.groundstate_itp(2e-2,1e-6);
% Calculate dynamics
dt = 0.01; % time step
totaliter = ceil(totaltime/(dt*100)); % total iteration number is defined by the total evolution time
figure; hold on;
Zres = zeros(1,15);
for i=1:15
    task.current_iter=0; % reset iteration counter for repeating calculations
    vv = i*(5e-8)/(rscale/tcoef); % barrier velocity
    tic
    task.solve_split(dt,100,totaliter);
    toc
    plot(task.history.Z); drawnow;
    Zres(i) = gather(task.history.Z(end));
end
figure; plot((1:15)*(5e-8),Zres,'+');