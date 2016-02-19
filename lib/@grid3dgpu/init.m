function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 3
  %  treat input as grid arrays
  x = gpuArray(varargin{ 1 });
  y = gpuArray(varargin{ 2 });
  z = gpuArray(varargin{ 3 });
else
  %  construct arrays
  x = gpuArray(linspace( -varargin{ 1 }, varargin{ 1 }, varargin{ 2 } ));
  y = gpuArray(linspace( -varargin{ 3 }, varargin{ 3 }, varargin{ 4 } ));
  z = gpuArray(linspace( -varargin{ 5 }, varargin{ 5 }, varargin{ 6 } ));
end

%  number of grid points
[ nx, ny, nz ] = deal( numel( x ), numel( y ), numel( z ) );
[ obj.nx, obj.ny, obj.nz ] = deal( nx, ny, nz );
%  total number of grid points
obj.n = nx * ny * nz;
%  spatial increment
dx = x( 2 ) - x( 1 );
dy = y( 2 ) - y( 1 );
dz = z( 2 ) - z( 1 );

%  save grid to object
[ obj.x, obj.y, obj.z ] = deal( x, y, z );

%%  wavevectors and operators in reciprocal space
kx = gpuArray([ (0:nx/2) -(nx/2-1:-1:1)]*pi/x(end));
ky = gpuArray([ (0:ny/2) -(ny/2-1:-1:1)]*pi/y(end));
kz = gpuArray([ (0:nz/2) -(nz/2-1:-1:1)]*pi/z(end));

[ kx, ky, kz ] = meshgrid( kx, ky, kz );
obj.kk = (kx.^2+ky.^2+kz.^2)/2;

obj.kx = gpuArray(2 * pi * ( 0 : nx - 1 ) / nx);
obj.ky = gpuArray(2 * pi * ( 0 : ny - 1 ) / ny);
obj.kz = gpuArray(2 * pi * ( 0 : nz - 1 ) / nz);
%  make mesh
[ kx, ky, kz ] = meshgrid( obj.kx, obj.ky, obj.kz );
                     
%  Laplace operator in Fourier space
obj.ilap = - 2 * ( 1 - cos( kx ) ) / dx ^ 2 -  ...
             2 * ( 1 - cos( ky ) ) / dy ^ 2 -  ...
             2 * ( 1 - cos( kz ) ) / dz ^ 2;
%  derivatives in Fourier space
obj.igradx = 1i * sin( kx ) / dx;
obj.igrady = 1i * sin( ky ) / dy;
obj.igradz = 1i * sin( kz ) / dz;

%%  meshgrid coordinates
%  mesh grid data 
[ x, y, z ] = meshgrid( x, y, z );
[x2, y2] = meshgrid(obj.x, obj.y);
%  save meshgrid structure
obj.mesh = struct( 'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2 );

obj.weight = dx*dy*dz;
