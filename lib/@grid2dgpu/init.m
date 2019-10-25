function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 2
  %  treat input as grid arrays
  x = gpuArray(varargin{ 1 });
  y = gpuArray(varargin{ 2 });
else
  %  construct arrays
  x = gpuArray(linspace( -varargin{ 1 }, varargin{ 1 }, varargin{ 2 } ));
  y = gpuArray(linspace( -varargin{ 3 }, varargin{ 3 }, varargin{ 4 } ));
end

%  number of grid points
[ nx, ny] = deal( numel( x ), numel( y ));
[ obj.nx, obj.ny ] = deal( nx, ny );
%  spatial increment
dx = x( 2 ) - x( 1 );
dy = y( 2 ) - y( 1 );

%  save grid to object
[ obj.x, obj.y ] = deal( x, y );

%%  wavevectors and operators in reciprocal space
kx = gpuArray([ (0:nx/2) -(nx/2-1:-1:1)]*pi/x(end));
ky = gpuArray([ (0:ny/2) -(ny/2-1:-1:1)]*pi/y(end));

[ obj.kx, obj.ky ] = meshgrid( kx, ky );
obj.kk = (obj.kx.^2+obj.ky.^2)/2;

%%  meshgrid coordinates
%  mesh grid data 
[ x, y ] = meshgrid( x, y );
%  save meshgrid structure
obj.mesh = struct( 'x', x, 'y', y, 'z', 0, 'x2', x, 'y2', y );

obj.weight = dx*dy;
