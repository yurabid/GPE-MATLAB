function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 1
  %  treat input as grid arrays
  x = (varargin{ 1 });
else
  %  construct arrays
  x = (linspace( -varargin{ 1 }, varargin{ 1 }, varargin{ 2 } ));
end

%  number of grid points
nx = numel( x );
obj.nx = nx;
%  spatial increment
dx = x( 2 ) - x( 1 );

%  save grid to object
obj.x = x;

%%  wavevectors and operators in reciprocal space
obj.kx = ([ (0:nx/2) -(nx/2-1:-1:1)]*pi/x(end));

obj.kk = (obj.kx.^2)/2;

%  save meshgrid structure
obj.mesh = struct( 'x', x, 'y', 0, 'z', 0, 'x2', x, 'y2', 0 );

obj.weight = dx;
