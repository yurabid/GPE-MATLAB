classdef grid2dgpu < handle
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    n           %  NOT USED - number of positions
    nx          %  number of x positions
    ny          %  number of y positions
    x           %  positions of grid along x
    y           %  positions of grid along y
    weight      %  integration weight
    kk          %  Laplace operator in Fourier space
    kx          %  NOT USED - wavenumbers along x
    ky          %  NOT USED - wavenumbers along y
    igradx      %  NOT USED - Fourier transform of derivative operator along x
    igrady      %  NOT USED - Fourier transform of derivative operator along y
    ilap        %  NOT USED - Fourier transform of Laplace operator
    mesh        %  meshgrid coordinates for x, y (also z=0, x2=x, y2=y are included for compatibility)
  end
  
%%  Methods
  methods
    
    function obj = grid2dgpu( varargin )
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
      obj = init( obj, varargin{ : } );
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid3d :' );
      disp( struct( 'x', gather([ min( obj.x ), max( obj.x ), length( obj.x ) ]),  ... 
                    'y', gather([ min( obj.y ), max( obj.y ), length( obj.y ) ]) ) );
    end
    
    function [] = imagesc( obj, fun, varargin )
      %  Plot function on grid.
      imagesc( gather(obj.x), gather(obj.y), gather(fun), varargin{ : } );
    end
    
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
