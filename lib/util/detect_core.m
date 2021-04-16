function [resp,resm] = detect_core(phi,rx,ry,refine)
% Detect vortex cores on a 2D slice of the wave function
% IN:
%   - phi: 2D array of complex numbers (slice of the wave function)
%   - rx, ry: grid of coordinate ponts in x and y 
%   - refine: (optional) refine vortex positions to sub-grid precision (can be slow with many vortices) 
% OUT:
%   - resp, resm: 2D arrays, each row  represents one detected vortex and contain [x,y] of the core. 
%                 resp (resm) contains only positive (negative) charged vortices
if(nargin==3)
    refine=0;
end
[n, m] = size(phi);
dx = rx(1,end/2+1)-rx(1,end/2);
dy = ry(end/2+1,1)-ry(end/2,1);
an = angle(phi);
rep = real(phi);
imp = imag(phi);
ddphi = abs(del2(abs(phi).^2,dx,dy));
idx = [0 0; 1 0; 1 1; 0 1];
nshift = size(idx,1);
angs = zeros(n,m,nshift,'like',rx);
for i=1:nshift
	angs(:,:,i) = circshift(an,idx(i,:));
end
dif = angs - circshift(angs,[0 0 1]);
res1 = (sum(dif>pi,3)-sum(dif<-pi,3)).*(ddphi>0.01*max(ddphi(:)));

res = res1>0;
if(refine>0)
    [ix,iy]=find(res);
    resp = zeros(numel(ix),2);
    for ii=1:numel(ix)
        resp(ii,:) = refine_core(max(ix(ii),2),max(iy(ii),2));
    end
else
    resp = [rx(res)-dx/2 ry(res)-dy/2];
end
res = res1<0;
if(refine>0)
    [ix,iy]=find(res);
    resm = zeros(numel(ix),2);
    for ii=1:numel(ix)
        resm(ii,:) = refine_core(max(ix(ii),2),max(iy(ii),2));
    end
else
    resm = [rx(res)-dx/2 ry(res)-dx/2];
end

    function res = refine_core(ix,iy)
    % Sub-grid-precision vortex detection    
    % Algorithm based on C.L. Phillips et.al. Phys.Rev.E 91, 023311 (2015)
        A = -[rep(ix-1,iy)-rep(ix,iy),rep(ix,iy);imp(ix-1,iy)-imp(ix,iy),imp(ix,iy)];
        B = [rep(ix,iy)-rep(ix-1,iy)-rep(ix,iy-1)+rep(ix-1,iy-1),rep(ix,iy-1)-rep(ix,iy);...
            imp(ix,iy)-imp(ix-1,iy)-imp(ix,iy-1)+imp(ix-1,iy-1),imp(ix,iy-1)-imp(ix,iy)];

        [v,d] = eig(A,B);
        res = [0,0];
        v(1,:)=v(1,:)./v(2,:);
        for i=1:2
            if((d(i,i)>0)&&(d(i,i)<1) && (v(1,i)>0) && (v(1,i)<1))
                res= [rx(ix,iy)-dy*d(i,i),ry(ix,iy)-dx*v(1,i)];
                break
            end
        end
    end
end