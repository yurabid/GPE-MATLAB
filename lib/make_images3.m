
if(exist('L','var') ~= 1)
    disp('reinitializing config parameters');
    config
end
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
dt = ddt*tcoef*niter_inner;
load('params.mat'); % read LL and maxx
lng = numel(find(abs(LL)));
slice = zeros(128,128);
for j=1:lng
    load(sprintf('snapshots/slicexz_%05d.mat',j));
    %load(sprintf('snapshots/slicexz_%05d.mat',j));
    %figure4(r, r, slice, slicez.', (1:lng)*dt, j, LL(1:lng), 0, maxx);
    figure3(r, rz, slicez', (1:lng)*dt, j, LL(1:lng), 0, maxx);
end