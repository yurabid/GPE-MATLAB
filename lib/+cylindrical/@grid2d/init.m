function  obj = init( obj, varargin )
%  INIT - Initialize 3D grid.

%%  set up real space grid from input
if length( varargin ) == 2
  %  treat input as grid arrays
  obj.r = (varargin{1});
  rmin = obj.r(1);
  rmax = obj.r(end);
  obj.phi = varargin{2};
  obj.dphi = obj.phi(2) - obj.phi(1);
  [obj.nr, obj.nphi] = deal(numel(obj.r), numel(obj.phi));
else
  %  construct arrays
  rmin = varargin{1};
  rmax = varargin{2};
  obj.nr = varargin{3};
  obj.nphi = varargin{4};
  obj.r = (linspace( rmin, rmax, obj.nr ));
  phi = (linspace( -pi, pi, obj.nphi + 1 ));
  obj.dphi = phi(2) - phi(1);
  obj.phi = phi(1:end-1) + obj.dphi/2;
end
obj.dr = (obj.r(2)-obj.r(1));
%%  wavevectors and operators in reciprocal space
kr = ([ (0:obj.nr/2) -(obj.nr/2-1:-1:1)]*2*pi/(rmax-rmin+obj.dr));
kphi = ([ (0:obj.nphi/2) -(obj.nphi/2-1:-1:1)]);

[ obj.kr, obj.kphi ] = meshgrid( kr, kphi );

%%  meshgrid coordinates
%  mesh grid data 
[ r, phi ] = meshgrid( obj.r, obj.phi );
obj.x = obj.r;
%  save meshgrid structure
obj.mesh = struct( 'r', r, 'phi', phi, 'z', 0, 'x2', r.*cos(phi), 'y2', r.*sin(phi), 'x', r, 'y', 0 );

obj.weight = obj.dr*obj.dphi;
obj.size = size(obj.mesh.x);
end
