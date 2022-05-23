function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 3
  %  treat input as grid arrays
  r = (varargin{ 1 });
  phi = (varargin{ 2 });
  z = (varargin{ 3 });
  hw = (r(end)-r(1))/2;
  [ nr, nphi, nz] = deal( numel( r ), numel( phi ), numel(z));
else
  %  construct arrays
  r0 = varargin{ 1 };
  hw = varargin{ 2 };
  r = (linspace( r0-hw, r0+hw, varargin{ 3 } ));
  phi = (linspace( -pi, pi, varargin{ 4 } + 1 ));
  dphi = phi( 2 ) - phi( 1 );
  phi = phi(1:end-1) + dphi/2;
  z = (linspace( -varargin{ 5 }, varargin{ 5 }, varargin{ 6 } ));
  nr = varargin{ 3 };
  nphi = varargin{ 4 };
  nz = varargin{ 6 };
end

[ obj.nr, obj.nphi, obj.nz ] = deal( nr, nphi, nz );
%  spatial increment
dr = r( 2 ) - r( 1 );
dphi = phi( 2 ) - phi( 1 );
dz = z( 2 ) - z( 1 );

%  save grid to object
[ obj.r, obj.phi, obj.z ] = deal( r, phi, z );

%%  wavevectors and operators in reciprocal space
kr = ([ (0:nr/2) -(nr/2-1:-1:1)]*pi/(hw+dr/2));
kphi = ([ (0:nphi/2) -(nphi/2-1:-1:1)]);
kz = ([ (0:nz/2) -(nz/2-1:-1:1)]*2*pi/(z(end)-z(1)+dz));

[ obj.kr, obj.kphi, obj.kz ] = meshgrid( kr, kphi, kz );
% obj.kk = (obj.kx.^2+obj.ky.^2)/2;

%%  meshgrid coordinates
%  mesh grid data 
[ r2, phi2] = meshgrid( r, phi);
[ r, phi, z ] = meshgrid( r, phi, z );
%  save meshgrid structure
obj.mesh = struct( 'r', r, 'phi', phi, 'z', z, 'x2', r2.*cos(phi2), 'y2', r2.*sin(phi2) );

obj.weight = dr*dphi*dz;
