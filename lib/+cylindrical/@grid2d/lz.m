function ret = lz(obj,phi)
	ret = obj.ifftphi(obj.kphi.*obj.fftphi(phi));
end
