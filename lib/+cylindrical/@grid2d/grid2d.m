classdef grid2d < grid2d
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
	nr
	nphi
	r
	phi
	kr
	kphi
    dr
    dphi
  end
  
%%  Methods
  methods
    
    function obj = grid2d( varargin )
      %  Initialize 2D polar grid.
      %
      %  Usage :
      %    obj = cylindrical.grid2d( r, phi)
      %    obj = cylindrical.grid2d( rmin, rmax, nr, nphi)
      %  Input
      %    r,    phi    :  arrays with grid points
      %    rmin, rmax   :  radial grid spans between rmin and rmax
      %    nr,   nphi   :  number of points in each dimension
      %  initialization
      if nargin>0
        obj.ndims = 2;
        obj = init( obj, varargin{ : } );
      end
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid2d :' );
      disp( struct( 'r', ([ min( obj.r ), max( obj.r ), length( obj.r ) ]),  ... 
                    'phi', ([ length( obj.y ) ]) ) );
    end
  
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
