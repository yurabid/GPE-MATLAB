function ifun = fftr( obj, fun )
%  FFT - Fast Fourier transform on grid.
%
%  Usage for obj = grid3d :
%    ifun = fft( obj, fun )
%  Input
%    fun    :  function values in real space
%  Output
%    ifun   :  function values in reciprocal space

ifun = fft(fft(fun, obj.nr, 2), obj.nz, 3);
