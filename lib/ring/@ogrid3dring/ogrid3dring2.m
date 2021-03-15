classdef ogrid3dring2 < handle
  %  Three-dimensional grid for solution of non-linear Schroedinger equation.
  
%%  Properties
  properties
    ndims=3
	omr         % harmonic trap frequency in r direction
	r0          % ring radius
	omz         % harmonic trap frequency in z direction
	ecut        % highest energy of the C region
	nstates     % total number of basis states
    npoints     % total number of quadrature grid points
    nr          %  number of r positions
    nphi        %  number of phi positions
    nz          %  number of z positions
    nsr         %  number of osc states in r
    nsphi       %  number of pw  states in phi
    nsz         %  number of osc states in z
    x
    r           %  positions of grid along r
    phi         %  positions of grid along phi
    z           %  positions of grid along z
    wr          %  Gauss-Hermite integration weights along r
    wphi        %  integration weights along phi
    wz          %  Gauss-Hermite integration weights along z
    transr      %  transformation matrix along r
    transphi    %  transformation matrix along phi
    transz      %  transformation matrix along z
    rmat
    drmat
    r2mat
    rdrmat
    rdphimat
    kphi
    weight      %  integration weight
    kk          %  Laplace operator in Fourier space
    mesh        %  meshgrid coordinates for x, y, z and x2, y2 (2D mesh in XY plane)
    wtot        %  [nx*ny*nz] array of weights for coordinate-space integration
    etot        %  [nsx*nsy*nsz] array energies of basis states
    etotphi     %  [nsx*nsy*nsz] array energies of basis states
    mask        %  [nsx*nsy*nsz] boolean array to cut states above ecut
  end
  
%%  Methods
  methods
    
    function obj = ogrid3dring2( omr, r0, omz, ecut, grid_factor )
      %  Initialize 3D grid.
      %
      %  Usage :
      %    obj = ogrid3d( omr, r0, omz, ecut )
      %  Input
      %  initialization
		obj.omr = omr;
		obj.r0 = r0;
		obj.omz = omz;
		obj.ecut = ecut;
		obj.nsr = gather(ceil((ecut+0.5)/omr));
        pindmax = gather(ceil(sqrt(2*ecut)*r0));
		obj.nsphi = 2*pindmax+1;
%         obj.nsphi = 1;
		obj.nsz = gather(ceil((ecut+0.5)/omz));
        
        obj.nr = 2*(obj.nsr-1);
        obj.nphi = 1*(obj.nsphi);
        obj.nz = 2*(obj.nsz-1);
        
        if(nargin <= 4)
            grid_factor = 1;
        end
        
        [obj.r, obj.wr] = obj.gauss_hermite_wrap(obj.nr,grid_factor*omr);
        obj.r = obj.r + r0;
        obj.x = obj.r;
        obj.wr = obj.wr.*obj.r.*(obj.r>0);
%         [obj.y, obj.wy] = obj.gauss_hermite_wrap(obj.ny,grid_factor*omy);
        obj.phi = linspace(0,2*pi,obj.nphi+1);
        obj.phi = obj.phi(1:end-1);
        obj.wphi = (obj.phi(2)-obj.phi(1));
        [obj.z, obj.wz] = obj.gauss_hermite_wrap(obj.nz,grid_factor*omz);
        
%         obj.wx = obj.wx.*(obj.wx>1e-30);
%         obj.wy = obj.wy.*(obj.wy>1e-30);
%         obj.wz = obj.wz.*(obj.wz>1e-30);
        
        obj.transr = obj.trans2r(obj.r-r0,obj.nsr,obj.omr);%/sqrt(r0);
        obj.transphi = obj.trans2plane(obj.phi,pindmax);
        obj.transz = obj.trans2(obj.z,obj.nsz,obj.omz);
        obj.rmat = (spdiags([sqrt((0:obj.nsr-1)'/2),(1:obj.nsr)'*0+r0,sqrt((1:obj.nsr)'/2)],[1,0,-1],obj.nsr,obj.nsr))^(-1);
        obj.drmat = (spdiags([-sqrt((0:obj.nsr-1)'/2),sqrt((1:obj.nsr)'/2)],[1,-1],obj.nsr,obj.nsr));
        obj.r2mat = obj.rmat^2;
        obj.rdrmat = obj.rmat*obj.drmat;
		obj.nstates = obj.nsr*obj.nsphi*obj.nsz;
        obj.npoints = obj.nr*obj.nphi*obj.nz;

%         [y,x,z] = meshgrid(obj.y,obj.x,obj.z);
%         [y2, x2] = meshgrid(obj.y, obj.x);
%         [wy,wx,wz] = meshgrid(obj.wy,obj.wx,obj.wz);
%         obj.mesh = struct( 'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2, 'wx', wx, 'wy', wy, 'wz', wz);
        [phi,r,z] = meshgrid(obj.phi,obj.r,obj.z);
        x = r.*cos(phi);
        y = r.*sin(phi);
        [phi2,r2] = meshgrid(obj.phi,obj.r);
        x2 = r2.*cos(phi2);
        y2 = r2.*sin(phi2);
        [wphi,wr,wz] = meshgrid(ones(1,obj.nphi)*obj.wphi,obj.wr,obj.wz);

        obj.mesh = struct( 'r',r,'phi',phi,'x', x, 'y', y, 'z', z, 'x2', x2, 'y2', y2, 'wr', wr, 'wphi', wphi, 'wz', wz);

        obj.wtot = wr.*wphi.*wz; %.*exp(grid_factor*(omx*x.^2+omy*y.^2+omz*z.^2));
        [i,j,k] = ind2sub([obj.nsr, obj.nsphi, obj.nsz],(1:obj.nstates));
        obj.etot = (obj.to3d(omr*(i-0.5) + omz*(k-0.5)));% + (obj.to3d( (j-pindmax-1).^2/2 ));
        obj.kphi = (obj.to3d( (j-pindmax-1) ));
        obj.etotphi = (obj.to3d( (j-pindmax-1).^2/2 ));
        obj.mask = obj.etot+(obj.etotphi-1/8)/r0^2 <= ecut+0.5;
        [i,j,k] = ind2sub([obj.nr, obj.nsphi, obj.nsz],(1:obj.nr*obj.nsphi*obj.nsz));
%         obj.rdphimat = reshape(((j-pindmax-1).^2-1/4)/2./r0^2,[obj.nr, obj.nsphi, obj.nsz]);
        obj.rdphimat = reshape(((j-pindmax-1).^2-1/4)/2./obj.r(i).^2,[obj.nr, obj.nsphi, obj.nsz]);
    end
  
    function disp( obj )
      %  Command window display.
      disp( 'grid3d :' );
      disp( struct( 'r', ([ obj.nsr ]),  ... 
                    'phi', ([ obj.nsphi ]),  ...
                    'z', ([ obj.nsz ]) ) );
    end
    function [x,w] = gauss_hermite_wrap(obj,n,alpha)
        [x,w] = gauss_hermite(n);
        if(isa(obj.r0,'gpuArray'))
            x = gpuArray(x./sqrt(alpha));
            w = gpuArray(w./sqrt(alpha));
        else
            x = (x./sqrt(alpha));
            w = (w./sqrt(alpha));
        end
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
        if(isa(obj.r0,'gpuArray'))
            res = gpuArray(zeros(length(x),n));
            res(:,1) = gpuArray(obj.osc_state(0,gather(x),gather(om)));
            res(:,2) = gpuArray(obj.osc_state(1,gather(x),gather(om)));
        else
            res = zeros(length(x),n);
            res(:,1) = obj.osc_state(0,x,om);
            res(:,2) = obj.osc_state(1,x,om);
        end
        for i=2:n-1
            res(:,i+1) = sqrt(2*om./i).*x(:).*res(:,i) - sqrt((i-1)/i).*res(:,i-1);
        end
    end  
    function res = trans2r(obj,x,n,om)
        
        if(isa(obj.r0,'gpuArray'))
            res = gpuArray(zeros(length(x),n));
            res(:,1) = gpuArray(obj.osc_state(0,gather(x),gather(om)))./sqrt(x(:)+obj.r0);
            res(:,2) = gpuArray(obj.osc_state(1,gather(x),gather(om)))./sqrt(x(:)+obj.r0);
        else
            res = zeros(length(x),n);
            res(:,1) = obj.osc_state(0,x,om)./sqrt(x(:)+obj.r0);
            res(:,2) = obj.osc_state(1,x,om)./sqrt(x(:)+obj.r0);
        end
        for i=2:n-1
            res(:,i+1) = sqrt(2*om./i).*x(:).*res(:,i) - sqrt((i-1)/i).*res(:,i-1);
        end
    end      
    function res = trans2plane(obj,x,n)
        res = zeros(length(x),2*n+1);
        for i=-n:n
            res(:,i+n+1) = exp(1i*i*x)/sqrt(2*pi);
        end
    end      
    function res = trans1mom(obj,grid,n,om)
        nosc = [0:n-1];
        nn = repmat(nosc,length(grid),1);
        res = obj.osc_state(nosc,grid,om).*(-1i).^nn;
    end  
    function res = rdr(obj,phi)
        res = reshape(permute(phi,[3 2 1]),[obj.nsphi*obj.nsz,obj.nsr])*obj.rdrmat.';
        res = permute(reshape(res,[obj.nsz,obj.nsphi,obj.nsr]),[3 2 1]).*obj.mask;
    end
    function res = r2(obj,phi)
        res = reshape(permute(phi,[3 2 1]),[obj.nsphi*obj.nsz,obj.nsr])*obj.r2mat.';
        res = permute(reshape(res,[obj.nsz,obj.nsphi,obj.nsr]),[3 2 1]).*obj.mask;
    end    
    function res = ifftr(obj,phi)
        res = reshape(permute(phi,[3 2 1]),[obj.nsphi*obj.nsz,obj.nsr])*obj.transr.';
        res = permute(reshape(res,[obj.nsz,obj.nsphi,obj.nr]),[3 2 1]);
    end
    function res = fftr(obj,phi)
        res = reshape(permute(phi,[3 2 1]),[obj.nsphi*obj.nsz,obj.nr])*obj.transr;
        res = permute(reshape(res,[obj.nsz,obj.nsphi,obj.nsr]),[3 2 1]).*obj.mask;
    end   
    function res = ifftzphi(obj,phi)
        res = reshape(permute(phi,[1 3 2]),[obj.nr*obj.nsz,obj.nsphi])*obj.transphi';       
        res = reshape(permute(reshape(res,[obj.nr,obj.nsz,obj.nphi]),[1 3 2]),[obj.nr*obj.nphi,obj.nsz])*obj.transz.';
        res = reshape(res,[obj.nr,obj.nphi,obj.nz]);
    end
%     function res = fftzphi(obj,phi)
%         res = reshape(permute(phi.*obj.wtot,[1 3 2]),[obj.nr*obj.nz,obj.nphi])*obj.transphi;        
%         res = reshape(permute(reshape(res,[obj.nr,obj.nz,obj.nsphi]),[1 3 2]),[obj.nr*obj.nsphi,obj.nz])*obj.transz;
%         res = reshape(res,[obj.nr,obj.nsphi,obj.nsz]);
%     end 
%     function res = ifftzphi(obj,phi)
%         res = reshape(phi,[obj.nr*obj.nsphi,obj.nsz])*obj.transz.';
%         res = reshape(permute(reshape(res,[obj.nr,obj.nsphi,obj.nz]),[1 3 2]),[obj.nr*obj.nz,obj.nsphi])*obj.transphi';       
%         res = reshape(permute(res,[1 3 2]),[obj.nr,obj.nphi,obj.nz]);
%     end
    function res = fftzphi(obj,phi)
        res = reshape(phi.*obj.wtot,[obj.nr*obj.nphi,obj.nz])*obj.transz;        
        res = reshape(permute(reshape(res,[obj.nr,obj.nphi,obj.nsz]),[1 3 2]),[obj.nr*obj.nsz,obj.nphi])*obj.transphi;
        res = permute(reshape(res,[obj.nr,obj.nsz,obj.nsphi]),[1 3 2]);
    end         
    function res = grid2sp(obj,phi)
        res = obj.grid2sp_inner(phi.*obj.wtot);
    end
    function res = applyh0(obj,phi)
        res = obj.etot.*phi + obj.fftr(obj.rdphimat.*obj.ifftr(phi));
    end
    function res = grid2sp_inner(obj,phi)
        res = reshape(permute(phi,[2 3 1]),[obj.nphi*obj.nz,obj.nr])*obj.transr;
        res = reshape(permute(reshape(res,[obj.nphi,obj.nz,obj.nsr]),[3 2 1]),[obj.nsr*obj.nz,obj.nphi])*obj.transphi;        
        res = reshape(permute(reshape(res,[obj.nsr,obj.nz,obj.nsphi]),[1 3 2]),[obj.nsr*obj.nsphi,obj.nz])*obj.transz;
        res = reshape(res,[obj.nsr,obj.nsphi,obj.nsz]).*obj.mask;
    end

    function res = sp2grid(obj,phi)
        res = reshape(permute(phi,[2 3 1]),[obj.nsphi*obj.nsz,obj.nsr])*obj.transr.';
        res = reshape(permute(reshape(res,[obj.nsphi,obj.nsz,obj.nr]),[3 2 1]),[obj.nr*obj.nsz,obj.nsphi])*obj.transphi';       
        res = reshape(permute(reshape(res,[obj.nr,obj.nsz,obj.nphi]),[1 3 2]),[obj.nr*obj.nphi,obj.nsz])*obj.transz.';
        res = reshape(res,[obj.nr,obj.nphi,obj.nz]);
    end
    function res = sp2grid_arb(obj,f,r,phi,z)
        dim = [length(r) length(phi) length(z)];
        res = zeros(dim);
        pindmax = round((obj.nsphi-1)/2);
        j=pindmax+1;
        for k = 1:obj.nsz
            for j = 1:obj.nsphi
                for i = 1:obj.nsr
                    res = res + f(i,j,k).*obj.osc_state(i-1,r-obj.r0,obj.omr)./sqrt(r).'.*exp(1i*(j-pindmax-1)*phi)/sqrt(2*pi).*obj.osc_state(k-1,z,obj.omz);
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
            res = reshape(phi,[obj.nsr, obj.nsphi, obj.nsz]);
        else
            res = reshape(phi,[obj.nr, obj.nphi, obj.nz]);
        end
    end
    function res=normalize(obj,phi)
        res = phi./sqrt(sum(abs(obj.to1d(phi)).^2));
    end
    function res=lz(obj,f)
        res = obj.kphi.*f;
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
    function res=get2ddens(obj,phi)
        res = sum(phi,3);
    end
  end
  
  methods (Access = private)
    obj = init( obj, varargin );
  end
  
end
