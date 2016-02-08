% Main configuration of parameters
% USAGE: place this file together with TDPot.m and empty snapshots
% directory; run "itp_run" to get initial state, "splitstep_run" to
% simulate dynamics, "make_images3(ddt*tcoef*niter_inner)" to build images.
% Run "clear all" after each change to this file.

global tcoef Ub omegaBar rscale linit vv;
%   Box size
L = 25;
Lz = 25;
%   Grid size
N=128;
Nz=128;
%   Physical parameters

NN0 = 5.0e4; % number of particles 
hbar=1.054571800e-34;
mass=1.41922608e-25;
omegascale = 50*2*pi;
tcoef = 1/(omegascale); %time scale
rscale = sqrt(hbar/mass/omegascale);
g = 4*pi*58.19e-10/rscale;
gam = 0.00; % dissipation constant gamma
tau = 1000000 / tcoef; % decay constant
omega = 0.0; %rotation speed
omegaBar = 0.0;
linit=0;
vv = (5e-7)/(rscale/tcoef);

% Potential (time-independent part)
om = 0.5*(17.68/50)^2;
Ub=650/omegascale*2*pi; %*1.5;
% Vfun = @(X,Y,Z) gpuArray(om*Z.^2 + 0.5*Y.^2 + 0.5*(abs(X)-r0).^2);
Vfun = @(X,Y,Z) gpuArray(om*Z.^2 + om*Y.^2 + 0.5*X.^2 );
% Vfun = @(X,Y,Z) gpuArray(40.23275*(1-exp(-om/40.23275*Z.^2)) + 2.743142*(1-exp(-0.5/2.743142*(sqrt(X.^2+Y.^2)-r0).^2)) );

% ITP parameters
niter = 1000; % number of ITP itrerations
dt_itp = 0.02; % time step for ITP

% Dynamics parameters
start = 0;
ddt = 0.009973310011396/2; % time step for dynamics
niter_inner = 100; %number of internal iterations
niter_outer = 1000; %630; %number of external iterations
n_cn = 10; % numder of CN iterations for L calculation
saveSlices = [0 0 1]; % save wave function slices in [(yz) (xz) (xy)] planes
detectCores = [0 0 0]; % detect cores in [(yz) (xz) (xy)] planes
useTDPot = 1; % use time-dependent potential (must be provided is a separate TDPot.m)
usePostProcess = 0; % use post-processing function PostProcess.m after each iteration

% Modification of the GS before the dynamics simulation (imprint vortex etc.)
% s = 1;
% sr = 1;
% rl0 = 0.5;
% rr0 = 6;
% z0 = 0.1;
% xi=0.2;
% phi_mod = @(X,Y,Z) tanh(sqrt((X-rl0).^2+Y.^2)/xi).^s.*exp(-1i*s*atan2(X - rl0,Y)).* ... % imprint s-charged off-center vortex
%     tanh(sqrt((sqrt(RS)-rr0).^2+(Z-z0).^2)/xi).^sr.*exp(-1i*sr*atan2(sqrt(RS)-rr0,Z-z0)); % imprint vortex ring
