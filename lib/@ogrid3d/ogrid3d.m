classdef ogrid3d < handle
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
	omx
	omy
	omz
	ecut
	nstates
    npoints
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
    transmat
    weight      %  integration weight
    kk          %  Laplace operator in Fourier space
    mesh        %  meshgrid coordinates for x, y, z and x2, y2 (2D mesh in XY plane)
    wtot
    etot
    mask
  end
  
%%  Methods
  methods
    
    function obj = ogrid3d( omx, omy, omz, ecut )
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
        [obj.x, obj.wx] = obj.gauss_hermite_wrap(obj.nx,2);
        [obj.y, obj.wy] = obj.gauss_hermite_wrap(obj.ny,2);
        [obj.z, obj.wz] = obj.gauss_hermite_wrap(obj.nz,2);
        obj.transx = obj.trans1(obj.x,obj.nsx,obj.omx);
        obj.transy = obj.trans1(obj.y,obj.nsy,obj.omy);
        obj.transz = obj.trans1(obj.z,obj.nsz,obj.omz);
		obj.nstates = obj.nsx*obj.nsy*obj.nsz;
        obj.npoints = obj.nx*obj.ny*obj.nz;
%         obj.transmat = kron(obj.transz,kron(obj.transy,obj.transx));
        [y,x,z] = meshgrid(obj.y,obj.x,obj.z);
        [y2, x2] = meshgrid(obj.y, obj.x);
        obj.mesh = struct( 'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2);
        [wy,wx,wz] = meshgrid(obj.wy,obj.wx,obj.wz);
%         obj.mesh.wx = wx;
%         obj.mesh.wy = wy;
%         obj.mesh.wz = wz;
        obj.wtot = (wx.*wy.*wz.*exp(2*(x.^2+y.^2+z.^2)));
        [i,j,k] = ind2sub([obj.nsx, obj.nsy, obj.nsz],(1:obj.nstates));
        obj.etot = (obj.to3d(omx*(i-0.5) + omy*(j-0.5) + omz*(k-0.5)));
        obj.mask = obj.etot <= ecut+0.5;
%         obj.etot = kron(kron((1:obj.nsx)-0.5,(1:obj.nsy)-0.5),(1:obj.nsz)-0.5)';
%         [x, wx] = obj.gauss_hermite_wrap(obj.nx);
%         [y, wy] = obj.gauss_hermite_wrap(obj.ny);
%         [z, wz] = obj.gauss_hermite_wrap(obj.nz);        
%         [x,y,z] = meshgrid(x,y,z);
% %         obj.mesh = struct( 'x', x, 'y', y, 'z', z);
%         [wx,wy,wz] = meshgrid(wx,wy,wz);
% %         obj.mesh.wx = wx;
% %         obj.mesh.wy = wy;
% %         obj.mesh.wz = wz;
%         obj.wtot1 = obj.to1d(wx.*wy.*wz);
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
    
    function res = trans1(obj,grid,n,om)
        nosc = [0:n-1];
        res = obj.osc_state(nosc,grid,om);
    end
    function res = trans1mom(obj,grid,n,om)
        nosc = [0:n-1];
        nn = repmat(nosc,length(grid),1);
        res = obj.osc_state(nosc,grid,om).*(-1i).^nn;
    end    
    function res = grid2spop(obj,phi)
        res = obj.grid2sp(phi.*obj.wtot);
    end
    function res = grid2sp(obj,phi)
%         res = (phi);
%         res1 = zeros(obj.nx,obj.nsy,obj.nz,'like',obj.x);
%         for i=1:obj.nz
%             res1(:,:,i) = phi(:,:,i)*obj.transy;
%         end
%         res1 = permute(res1,[2 1 3]);        
        res2 = zeros(obj.nsx,obj.nsy,obj.nz,'like',obj.x);
        for i=1:obj.nz
            res2(:,:,i) = obj.transx.'*phi(:,:,i)*obj.transy;
        end        
        res2 = permute(res2,[1 3 2]);
        res = zeros(obj.nsx,obj.nsz,obj.nsy,'like',obj.x);
        for i=1:obj.nsy
            res(:,:,i) = res2(:,:,i)*obj.transz;
        end        
        res = permute(res,[1 3 2]).*obj.mask;
%         res = obj.to1d(res);
    end
%     function res = grid2sp_old(obj,phi)
%         res = obj.transmat.'*obj.to1d(phi);
%     end
%     function res = grid2spop_old(obj,phi)
%         res = obj.transmat.'*(phi.*obj.wtot);
%     end
    function res = sp2grid(obj,phi)
%         res = (phi);
        res1 = zeros(obj.nx,obj.ny,obj.nsz,'like',obj.x);
        for i=1:obj.nsz
            res1(:,:,i) = obj.transx*phi(:,:,i)*obj.transy.';
        end
        res1 = permute(res1,[1 3 2]);        
%         res2 = zeros(obj.ny,obj.nx,obj.nsz,'like',obj.x);
%         for i=1:obj.nsz
%             res2(:,:,i) = res1(:,:,i)*obj.transx.';
%         end        
%         res2 = permute(res2,[2 3 1]);
        res = zeros(obj.nx,obj.nz,obj.ny,'like',obj.x);
        for i=1:obj.ny
            res(:,:,i) = res1(:,:,i)*obj.transz.';
        end        
        res = permute(res,[1 3 2]);
%         res = obj.to1d(res);        
%         res = obj.transmat*phi;
    end
    function res = sp2grid_arb(obj,phi,x,y,z)
        dim = [length(x) length(y) length(z)];
        res = zeros(dim);
        for k = 1:obj.nsz
            for j = 1:obj.nsy
                for i = 1:obj.nsx
%                     tmp = obj.osc_state(i-1,x,obj.omx)*obj.osc_state(j-1,y,obj.omy).';
%                     tmp2 = reshape(tmp,[],1);
%                     tmp3 = tmp2*obj.osc_state(k-1,z,obj.omz).';
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
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
