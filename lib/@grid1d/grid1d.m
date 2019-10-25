classdef grid1d < gridBase
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    nx          %  number of x positions
    x           %  positions of grid along x
    kx
%     weight      %  integration weight
%     kk          %  Laplace operator in Fourier space
%     mesh        %  meshgrid coordinates included for compatibility
  end
  
%%  Methods
  methods
    
    function obj = grid1d( varargin )
      %  Initialize 3D grid.
      %
      %  Usage :
      %    obj = grid1d( x )
      %    obj = grid1d( xmax, nx)
      %  Input
      %    x      :  array with grid points
      %    xmax   :  largest  grid points (smallest points are symmetric)
      %    nx     :  number of points 
      %  initialization
      obj.ndims = 1;
      obj = init( obj, varargin{ : } );
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid1d :' );
      disp( struct( 'x', ([ min( obj.x ), max( obj.x ), length( obj.x ) ]) ) );
    end
    
    function [] = imagesc( obj, fun, varargin )
      %  Plot function on grid.
      plot( (obj.x), (fun), varargin{ : } );
    end
    
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
