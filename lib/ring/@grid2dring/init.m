function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 2
  %  treat input as grid arrays
  r = (varargin{ 1 });
  phi = (varargin{ 2 });
  hw = (r(end)-r(1))/2;
  [ nr, nphi] = deal( numel( r ), numel( phi ));
else
  %  construct arrays
  r0 = varargin{ 1 };
  hw = varargin{ 2 };
  r = (linspace( r0-hw, r0+hw, varargin{ 3 } ));
  phi = (linspace( -pi, pi, varargin{ 4 } + 1 ));
  dphi = phi( 2 ) - phi( 1 );
  phi = phi(1:end-1) + dphi/2;
  nr = varargin{ 3 };
  nphi = varargin{ 4 };
end

[ obj.nr, obj.nphi ] = deal( nr, nphi );
%  spatial increment
dr = r( 2 ) - r( 1 );
%  dphi = phi( 2 ) - phi( 1 );

%  save grid to object
[ obj.r, obj.phi ] = deal( r, phi );

%%  wavevectors and operators in reciprocal space
kr = ([ (0:nr/2) -(nr/2-1:-1:1)]*pi/(hw+dr/2));
kphi = ([ (0:nphi/2) -(nphi/2-1:-1:1)]);

[ obj.kr, obj.kphi ] = meshgrid( kr, kphi );
% obj.kk = (obj.kx.^2+obj.ky.^2)/2;

%%  meshgrid coordinates
%  mesh grid data 
[ r, phi ] = meshgrid( r, phi );
%  save meshgrid structure
obj.mesh = struct( 'r', r, 'phi', phi, 'z', 0, 'x2', r.*cos(phi), 'y2', r.*sin(phi) );

obj.weight = dr*dphi;
