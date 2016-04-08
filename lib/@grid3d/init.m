function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 3
  %  treat input as grid arrays
  x = (varargin{ 1 });
  y = (varargin{ 2 });
  z = (varargin{ 3 });
else
  %  construct arrays
  x = (linspace( -varargin{ 1 }, varargin{ 1 }, varargin{ 2 } ));
  y = (linspace( -varargin{ 3 }, varargin{ 3 }, varargin{ 4 } ));
  z = (linspace( -varargin{ 5 }, varargin{ 5 }, varargin{ 6 } ));
end

%  number of grid points
[ nx, ny, nz ] = deal( numel( x ), numel( y ), numel( z ) );
[ obj.nx, obj.ny, obj.nz ] = deal( nx, ny, nz );
%  spatial increment
dx = x( 2 ) - x( 1 );
dy = y( 2 ) - y( 1 );
dz = z( 2 ) - z( 1 );

%  save grid to object
[ obj.x, obj.y, obj.z ] = deal( x, y, z );

%%  wavevectors and operators in reciprocal space
kx = ([ (0:nx/2) -(nx/2-1:-1:1)]*pi/x(end));
ky = ([ (0:ny/2) -(ny/2-1:-1:1)]*pi/y(end));
kz = ([ (0:nz/2) -(nz/2-1:-1:1)]*pi/z(end));

[ kx, ky, kz ] = meshgrid( kx, ky, kz );
obj.kk = (kx.^2+ky.^2+kz.^2)/2;

%%  meshgrid coordinates
%  mesh grid data 
[ x, y, z ] = meshgrid( x, y, z );
[x2, y2] = meshgrid(obj.x, obj.y);
%  save meshgrid structure
obj.mesh = struct( 'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2 );

obj.weight = dx*dy*dz;
