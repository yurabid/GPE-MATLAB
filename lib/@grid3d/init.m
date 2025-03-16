function  obj = init( obj, varargin )
%  Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 3
  %  treat input as grid arrays
  x = (varargin{ 1 });
  y = (varargin{ 2 });
  z = (varargin{ 3 });
  [ nx, ny, nz ] = deal( numel( x ), numel( y ), numel( z ) );
else
  %  construct arrays
  x = (linspace( -varargin{ 1 }, varargin{ 1 }, varargin{ 2 } ));
  y = (linspace( -varargin{ 3 }, varargin{ 3 }, varargin{ 4 } ));
  z = (linspace( -varargin{ 5 }, varargin{ 5 }, varargin{ 6 } ));
  nx = varargin{ 2 };
  ny = varargin{ 4 };  
  nz = varargin{ 6 };  
end

%  number of grid points

[ obj.nx, obj.ny, obj.nz ] = deal( nx, ny, nz );
%  spatial increment
obj.dx = x( 2 ) - x( 1 );
obj.dy = y( 2 ) - y( 1 );
obj.dz = z( 2 ) - z( 1 );

%  save grid to object
[ obj.x, obj.y, obj.z ] = deal( x, y, z );

%  wavevectors and operators in reciprocal space
kx = ([ (0:nx/2) -(nx/2-1:-1:1)]*2*pi/(x(end)-x(1)+obj.dx));
ky = ([ (0:ny/2) -(ny/2-1:-1:1)]*2*pi/(y(end)-y(1)+obj.dy));
kz = ([ (0:nz/2) -(nz/2-1:-1:1)]*2*pi/(z(end)-z(1)+obj.dz));

[ obj.kx, obj.ky, obj.kz ] = meshgrid( kx, ky, kz );
obj.kk = (obj.kx.^2+obj.ky.^2+obj.kz.^2)/2;

%  meshgrid coordinates
[ x, y, z ] = meshgrid( x, y, z );
[x2, y2] = meshgrid(obj.x, obj.y);
%  save meshgrid structure
obj.mesh = struct( 'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2 );
obj.weight = obj.dx*obj.dy*obj.dz;
