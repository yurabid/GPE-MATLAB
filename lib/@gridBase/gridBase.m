classdef gridBase < handle
  %  Base grid object. Not to be used directly
  
%%  Properties
  properties
   
    weight      %  integration weight
    kweight     %  integration weight in fourier space
    kk          %  Laplace operator in Fourier space
    mesh        %  meshgrid coordinates
    ndims       %  number of dimensions
  end
  
%%  Methods
  methods
    
    function obj = gridBase( varargin )

    end
  
    
  end
  
  
end
