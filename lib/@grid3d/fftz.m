function ifun = fftz( obj, fun )
%  FFT - Fast Fourier transform on the grid.
%
%  Usage for obj = grid3d :
%    ifun = fft( obj, fun )
%  Input
%    fun    :  function values in real space
%  Output
%    ifun   :  function values in reciprocal space

ifun = fft(fun, obj.nz, 3);
