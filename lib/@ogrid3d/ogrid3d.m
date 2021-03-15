classdef ogrid3d < handle
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
	omx         % harmonic trap frequency in x direction
	omy         % harmonic trap frequency in y direction
	omz         % harmonic trap frequency in z direction
	ecut        % highest energy of the C region
	nstates     % total number of basis states
    npoints     % total number of quadrature grid points
    nx          %  number of x positions
    ny          %  number of y positions
    nz          %  number of z positions
    nsx         %  number of osc states in x
    nsy         %  number of osc states in y
    nsz         %  number of osc states in z
    x           %  positions of grid along x
    y           %  positions of grid along y
    z           %  positions of grid along z
    wx          %  Gauss-Hermite integration weights along x
    wy          %  Gauss-Hermite integration weights along y
    wz          %  Gauss-Hermite integration weights along z
    transx      %  transformation matrix along x
    transy      %  transformation matrix along y
    transz      %  transformation matrix along z
    weight      %  integration weight
    kk          %  Laplace operator in Fourier space
    mesh        %  meshgrid coordinates for x, y, z and x2, y2 (2D mesh in XY plane)
    wtot        %  [nx*ny*nz] array of weights for coordinate-space integration
    etot        %  [nsx*nsy*nsz] array energies of basis states
    mask        %  [nsx*nsy*nsz] boolean array to cut states above ecut
  end
  
%%  Methods
  methods
    
    function obj = ogrid3d( omx, omy, omz, ecut, grid_factor )
      %  Initialize 3D grid.
      %
      %  Usage :
      %    obj = ogrid3d( omx, omy, omz, ecut )
      %  Input
      %  initialization
		obj.omx = omx;
		obj.omy = omy;
		obj.omz = omz;
		obj.ecut = ecut;
		obj.nsx = ceil((ecut+0.5)/omx);
		obj.nsy = ceil((ecut+0.5)/omy);
		obj.nsz = ceil((ecut+0.5)/omz);
        obj.nx = 2*(obj.nsx-1);
        obj.ny = 2*(obj.nsy-1);
        obj.nz = 2*(obj.nsz-1);
        
        if(nargin <= 4)
            grid_factor = 2;
        end
        
        [obj.x, obj.wx] = obj.gauss_hermite_wrap(obj.nx,grid_factor*omx);
        [obj.y, obj.wy] = obj.gauss_hermite_wrap(obj.ny,grid_factor*omy);
        [obj.z, obj.wz] = obj.gauss_hermite_wrap(obj.nz,grid_factor*omz);
        
%         obj.wx = obj.wx.*(obj.wx>1e-30);
%         obj.wy = obj.wy.*(obj.wy>1e-30);
%         obj.wz = obj.wz.*(obj.wz>1e-30);
        
        obj.transx = obj.trans2(obj.x,obj.nsx,obj.omx);
        obj.transy = obj.trans2(obj.y,obj.nsy,obj.omy);
        obj.transz = obj.trans2(obj.z,obj.nsz,obj.omz);
        
		obj.nstates = obj.nsx*obj.nsy*obj.nsz;
        obj.npoints = obj.nx*obj.ny*obj.nz;

        [y,x,z] = meshgrid(obj.y,obj.x,obj.z);
        [y2, x2] = meshgrid(obj.y, obj.x);
        [wy,wx,wz] = meshgrid(obj.wy,obj.wx,obj.wz);
        obj.mesh = struct( 'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2, 'wx', wx, 'wy', wy, 'wz', wz);

        obj.wtot = wx.*wy.*wz; %.*exp(grid_factor*(omx*x.^2+omy*y.^2+omz*z.^2));
        [i,j,k] = ind2sub([obj.nsx, obj.nsy, obj.nsz],(1:obj.nstates));
        obj.etot = (obj.to3d(omx*(i-0.5) + omy*(j-0.5) + omz*(k-0.5)));
        obj.mask = obj.etot <= ecut+0.5;
 
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid3d :' );
      disp( struct( 'x', ([ obj.nsx ]),  ... 
                    'y', ([ obj.nsy ]),  ...
                    'z', ([ obj.nsz ]) ) );
    end
    function [x,w] = gauss_hermite_wrap(obj,n,alpha)
        [x,w] = gauss_hermite(n);
        x = (x./sqrt(alpha));
        w = (w./sqrt(alpha));
    end
	function u=osc_state(obj,n,x,om)
        [~, xgrid] = meshgrid(n,x);
		u = hermite(n,x.*sqrt(om),'norm').*exp(-om*xgrid.^2*0.5)*om^0.25; %./sqrt(2^n*factorial(n)*sqrt(pi/om));
    end
    
	function u=osc_state2(obj,n,x,om)
        if(n<2)
            u = obj.osc_state(n,x,om);
        else
            u = sqrt(2*om./n).*x(:).*obj.osc_state2(n-1,x,om) - sqrt((n-1)/n).*obj.osc_state2(n-2,x,om);
        end
    end
    
    function res = trans1(obj,grid,n,om)
        nosc = [0:n-1];
        res = obj.osc_state(nosc,grid,om);
    end
    function res = trans2(obj,x,n,om)
%         nosc = [0:n-1];
%         res = obj.osc_state(nosc,grid,om);
        res = zeros(length(x),n);
        res(:,1) = obj.osc_state(0,x,om);
        res(:,2) = obj.osc_state(1,x,om);
        for i=2:n-1
            res(:,i+1) = sqrt(2*om./i).*x(:).*res(:,i) - sqrt((i-1)/i).*res(:,i-1);
        end
    end    
    function res = trans1mom(obj,grid,n,om)
        nosc = [0:n-1];
        nn = repmat(nosc,length(grid),1);
        res = obj.osc_state(nosc,grid,om).*(-1i).^nn;
    end
    function res = applyh0(obj,phi)
        res = obj.etot.*phi;
    end
    function res = grid2sp(obj,phi)
        res = obj.grid2sp_inner(phi.*obj.wtot);
    end
    function res = grid2sp_inner(obj,phi)
        res = reshape(permute(phi,[2 3 1]),[obj.ny*obj.nz,obj.nx])*obj.transx;
        res = reshape(permute(reshape(res,[obj.ny,obj.nz,obj.nsx]),[3 2 1]),[obj.nsx*obj.nz,obj.ny])*obj.transy;        
        res = reshape(permute(reshape(res,[obj.nsx,obj.nz,obj.nsy]),[1 3 2]),[obj.nsx*obj.nsy,obj.nz])*obj.transz;
        res = reshape(res,[obj.nsx,obj.nsy,obj.nsz]).*obj.mask;
    end

    function res = sp2grid(obj,phi)
        res = reshape(permute(phi,[2 3 1]),[obj.nsy*obj.nsz,obj.nsx])*obj.transx.';
        res = reshape(permute(reshape(res,[obj.nsy,obj.nsz,obj.nx]),[3 2 1]),[obj.nx*obj.nsz,obj.nsy])*obj.transy.';       
        res = reshape(permute(reshape(res,[obj.nx,obj.nsz,obj.ny]),[1 3 2]),[obj.nx*obj.ny,obj.nsz])*obj.transz.';
        res = reshape(res,[obj.nx,obj.ny,obj.nz]);
    end
    function res = sp2grid_arb(obj,phi,x,y,z)
        dim = [length(x) length(y) length(z)];
        res = zeros(dim);
        for k = 1:obj.nsz
            for j = 1:obj.nsy
                for i = 1:obj.nsx
                    res = res + phi(i,j,k)*reshape(reshape(obj.osc_state(i-1,x,obj.omx)*obj.osc_state(j-1,y,obj.omy).',[],1)*obj.osc_state(k-1,z,obj.omz).',dim);
                end
            end
        end
    end
    function res = to1d(obj,phi)
        res = phi(:);
    end
    function res = to3d(obj,phi)
        if(ndims(phi)==3)
            res = phi;
        elseif(length(phi)==obj.nstates)
            res = reshape(phi,[obj.nsx, obj.nsy, obj.nsz]);
        else
            res = reshape(phi,[obj.nx, obj.ny, obj.nz]);
        end
    end
    function res=normalize(obj,phi)
        res = phi./sqrt(sum(abs(obj.to1d(phi)).^2));
    end
    function res=inner(obj,f1,f2)
        res = sum(obj.to1d(conj(f1).*f2));
    end
    function res=integrate(obj,phi)
        res = sum(obj.to1d(phi));
    end
    function res=integrate_grid(obj,phi)
        res = sum(obj.to1d(phi.*obj.wtot));
    end
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
