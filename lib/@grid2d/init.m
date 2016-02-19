function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 2
  %  treat input as grid arrays
  x = (varargin{ 1 });
  y = (varargin{ 2 });
else
  %  construct arrays
  x = (linspace( -varargin{ 1 }, varargin{ 1 }, varargin{ 2 } ));
  y = (linspace( -varargin{ 3 }, varargin{ 3 }, varargin{ 4 } ));
end

%  number of grid points
[ nx, ny] = deal( numel( x ), numel( y ));
[ obj.nx, obj.ny ] = deal( nx, ny );
%  total number of grid points
obj.n = nx * ny;
%  spatial increment
dx = x( 2 ) - x( 1 );
dy = y( 2 ) - y( 1 );

%  save grid to object
[ obj.x, obj.y ] = deal( x, y );

%%  wavevectors and operators in reciprocal space
kx = ([ (0:nx/2) -(nx/2-1:-1:1)]*pi/x(end));
ky = ([ (0:ny/2) -(ny/2-1:-1:1)]*pi/y(end));

[ kx, ky ] = meshgrid( kx, ky );
obj.kk = (kx.^2+ky.^2)/2;

obj.kx = (2 * pi * ( 0 : nx - 1 ) / nx);
obj.ky = (2 * pi * ( 0 : ny - 1 ) / ny);
%  make mesh
[ kx, ky ] = meshgrid( obj.kx, obj.ky );
                     
%  Laplace operator in Fourier space
obj.ilap = - 2 * ( 1 - cos( kx ) ) / dx ^ 2 -  ...
             2 * ( 1 - cos( ky ) ) / dy ^ 2;
%  derivatives in Fourier space
obj.igradx = 1i * sin( kx ) / dx;
obj.igrady = 1i * sin( ky ) / dy;

%%  meshgrid coordinates
%  mesh grid data 
[ x, y ] = meshgrid( x, y );
%  save meshgrid structure
obj.mesh = struct( 'x', x, 'y', y, 'z', 0, 'x2', x, 'y2', y );

obj.weight = dx*dy;