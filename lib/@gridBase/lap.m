function ret = lap(obj,phi)
	ret = obj.ifft(obj.kk.*obj.fft(phi));
end
