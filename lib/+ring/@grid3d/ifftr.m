function fun = ifftr( obj, ifun )
%  FFT - Inverse fast Fourier transform on grid.
%
%  Usage for obj = grid3d :
%    fun = fft( obj, ifun )
%  Input
%    ifun   :  function values in reciprocal space
%  Output
%    fun    :  function values in real space

fun = ifft(ifft(ifun, obj.nz, 3), obj.nr, 2);
