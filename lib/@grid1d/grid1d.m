classdef grid1d < handle
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    nx          %  number of x positions
    x           %  positions of grid along x
    kx
    dx
    weight      %  integration weight
    kweight     %  integration weight in fourier space
    kk          %  Laplace operator in Fourier space
    mesh        %  meshgrid coordinates
    ndims       %  number of dimensions
    size        %  array of grid sizes in each dimension
    stencil1=[1,-8,0,8,-1]/12
    stencil2=[-1,16,-30,16,-1]/12
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
      if nargin>0
        obj.ndims = 1;
        obj = init( obj, varargin{ : } );
      end
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

    function df = deriv1(obj, A, dim, h)
        kernelSize = ones(1, ndims(A));
        kernelSize(dim) = 5;
        shapedKernel = reshape(obj.stencil1, kernelSize)/h;
        df = convn(A, shapedKernel, 'same');
    end

    function df = deriv2(obj, A, dim, h)
        kernelSize = ones(1, ndims(A));
        kernelSize(dim) = 5;
        shapedKernel = reshape(obj.stencil2, kernelSize)/h;
        df = convn(A, shapedKernel, 'same');
    end

  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
