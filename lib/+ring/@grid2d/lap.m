function ret = lap(obj,phi,omega)
if(nargin <= 2)
	omega = 0;
end
	ret = obj.ifftr((obj.kr.^2/2).*obj.fftr(phi)) + ...
		obj.ifftphi(((obj.kphi.^2-1/4)/2./obj.mesh.r.^2 + omega*obj.kphi).*obj.fftphi(phi));
end
