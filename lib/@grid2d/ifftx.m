function fun = ifftx( obj, ifun )
%  FFT - Inverse fast Fourier transform on grid.
%
%  Usage for obj = grid3d :
%    fun = fft( obj, ifun )
%  Input
%    ifun   :  function values in reciprocal space
%  Output
%    fun    :  function values in real space

fun = ifft(ifun, obj.nx, 2);
