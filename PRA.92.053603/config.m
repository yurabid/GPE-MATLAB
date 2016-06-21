% Main configuration of parameters
% USAGE: place this file together with TDPot.m and empty snapshots
% directory; run "itp_run" to get initial state, "splitstep_run" to
% simulate dynamics, "make_images3" to build images.
% Run "clear all" after each change to this file.

global tcoef;
%   Box size
L = 26;
Lz = 26;
%   Grid size
N=128;
Nz=128;
%   Physical parameters
g = 0.0276;
NN0 = 1.0e6; % number of particles 
tcoef = 1/(100*2*pi); %time scale
gam = 0.02; % dissipation constant gamma
tau = 100000 / tcoef; % decay constant
omega = 0.3; %rotation speed

% Potential (time-independent part)
r0 = 0;
om = 0.5*(5)^2;
vm = 34;
vz = 18;
alpha = 0.8;
beta = 0.2;
m = 1;
Vfun = @(X,Y,Z) gpuArray(0.5*Z.^2 * 9 + 0.5*(X.^2+Y.^2) + vm*beta^(2*m)*(X.^2+Y.^2).^m.*exp(-m*(beta^2*(X.^2+Y.^2) - 1)) + vz*exp(-alpha^2*Z.^2));
%Vfun = @(X,Y,Z) gpuArray(0.5*Z.^2 + 0.5*(X.^2+Y.^2)/4);

% ITP parameters
niter = 1000; % number of ITP itrerations
dt_itp = 0.004; % time step for ITP

% Dynamics parameters
start = 0;
ddt = 0.007; % time step for dynamics
niter_inner = 100; %number of internal iterations
niter_outer = 1000; %number of external iterations
n_cn = 10; % numder of CN iterations for L calculation
saveSlices = [0 1 0]; % save wave function slices in [(yz) (xz) (xy)] planes
detectCores = [1 1 1]; % detect cores in [(yz) (xz) (xy)] planes
useTDPot = 0; % use time-dependent potential (must be provided in a separate TDPot.m)

% Modification of the GS before the dynamics simulation (imprint vortex etc.)
    s = 1;
    sr = 1;
    rl0 = 0.5;
    rr0 = 6;
    z0 = 0.1;
    xi=0.1;
 phimod = @(X,Y,Z) tanh(sqrt((X-rl0).^2+Y.^2)/xi).^s.*exp(-1i*s*atan2(X - rl0,Y)).* ... % imprint s-charged off-center vortex
     tanh(sqrt((sqrt(X.^2+Y.^2)-rr0).^2+(Z-z0).^2)/xi).^sr.*exp(-1i*sr*atan2(sqrt(X.^2+Y.^2)-rr0,Z-z0)); % imprint vortex ring