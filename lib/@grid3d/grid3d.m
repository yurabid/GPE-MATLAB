classdef grid3d < grid2d
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    nz          %  number of z positions
    z           %  positions of grid along z
    kz          %  momentum grid along z      
  end
  
%%  Methods
  methods
    
    function obj = grid3d( varargin )
      %  Initialize 3D grid.
      %
      %  Usage :
      %    obj = grid3d( x, y, z )
      %    obj = grid3d( xmax, nx, ymax, ny, zmax, nz )
      %  Input
      %    x,    y,    z      :  arrays with grid points
      %    xmax, ymax, zmax   :  largest  grid points (smallest points are symmetric)
      %    nx,   ny,   nz     :  number of point in each dimension
      %  initialization
      if nargin>0
        obj.ndims = 3;
        obj = init( obj, varargin{ : } );
      end
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid3d :' );
      disp( struct( 'x', ([ min( obj.x ), max( obj.x ), length( obj.x ) ]),  ... 
                    'y', ([ min( obj.y ), max( obj.y ), length( obj.y ) ]),  ...
                    'z', ([ min( obj.z ), max( obj.z ), length( obj.z ) ]) ) );
    end
    
    function [] = imagesc( obj, fun, varargin )
      %  Plot function on grid.
      imagesc( (obj.x), (obj.y), (fun(:,:,obj.nz/2)), varargin{ : } );
    end
    
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
