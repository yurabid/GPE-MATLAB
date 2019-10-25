classdef grid3dgpu < gridBase
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    nx          %  number of x positions
    ny          %  number of y positions
    nz          %  number of z positions
    x           %  positions of grid along x
    y           %  positions of grid along y
    z           %  positions of grid along z
    kx          %  momentum grid along x
    ky          %  momentum grid along y
    kz          %  momentum grid along z    
%     weight      %  integration weight
%     kweight     %  integration weight in fourier space
%     kk          %  Laplace operator in Fourier space
%     mesh        %  meshgrid coordinates for x, y, z and x2, y2 (2D mesh in XY plane)
  end
  
%%  Methods
  methods
    
    function obj = grid3dgpu( varargin )
      %  Initialize 3D grid.
      %
      %  Usage :
      %    obj = grid3d( x, y, z )
      %    obj = grid3d( xmax, nx, ymax, ny, zmax, nz )
      %  Input
      %    x,    y,    z      :  arrays with grid points
      %    xmax, ymax, zmax   :  largest  grid points (smallest points are symmetric)
      %    nx,   ny,   nz     :  number of points in each dimension
      %  initialization
      obj.ndims = 3;
      obj = init( obj, varargin{ : } );
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid3d :' );
      disp( struct( 'x', gather([ min( obj.x ), max( obj.x ), length( obj.x ) ]),  ... 
                    'y', gather([ min( obj.y ), max( obj.y ), length( obj.y ) ]),  ...
                    'z', gather([ min( obj.z ), max( obj.z ), length( obj.z ) ]) ) );
    end
    
    function [] = imagesc( obj, fun, varargin )
      %  Plot function on grid.
      imagesc( gather(obj.x), gather(obj.y), gather(fun(:,:,obj.nz/2)), varargin{ : } );
    end
    
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
