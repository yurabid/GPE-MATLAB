classdef grid2dring < grid2d
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
	nr
	nphi
	r
	phi
	kr
	kphi
  end
  
%%  Methods
  methods
    
    function obj = grid2dring( varargin )
      %  Initialize 2D polar grid.
      %
      %  Usage :
      %    obj = grid3d( r, phi)
      %    obj = grid3d( r0, dr, nr, nphi)
      %  Input
      %    r,    phi    :  arrays with grid points
      %    r0, dr       :  radial grid spans between r0-dr and r0+dr
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
