classdef grid2d < grid1d
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    ny          %  number of y positions
    y           %  positions of grid along y
    ky
  end
  
%%  Methods
  methods
    
    function obj = grid2d( varargin )
      %  Initialize 3D grid.
      %
      %  Usage :
      %    obj = grid3d( x, y)
      %    obj = grid3d( xmax, nx, ymax, ny)
      %  Input
      %    x,    y      :  arrays with grid points
      %    xmax, ymax   :  largest  grid points (smallest points are symmetric)
      %    nx,   ny     :  number of point in each dimension
      %  initialization
      if nargin>0
        obj.ndims = 2;
        obj = init( obj, varargin{ : } );
      end
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid2d :' );
      disp( struct( 'x', ([ min( obj.x ), max( obj.x ), length( obj.x ) ]),  ... 
                    'y', ([ min( obj.y ), max( obj.y ), length( obj.y ) ]) ) );
    end
    
    function [] = imagesc( obj, fun, varargin )
      %  Plot function on grid.
      imagesc( (obj.x), (obj.y), (fun), varargin{ : } );
    end
    
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
