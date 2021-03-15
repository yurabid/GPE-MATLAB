function [frequencies, amplitude] = FSpectrum( signalpts, period )
%FSPECTRUM gives the single-sided Fourier-amplitudes and frequencies
%   Input :
%       signalpts   -   Discrete points of the signal; last point 
%                       is omitted when length(signalpts) is odd
%       period      -   time difference between two points in real time
%   Output :
%       amplitude   -   amplitude of the Fourier-spectrum
%       frequencies -   frequencies of the Fourier-spectrum in Hz (careful
%                       with 2*pi!)

dots = length(signalpts);
if mod(dots,2)==1
    dots=dots-1;
end

FTSignal = fft(signalpts);
amplitudedouble = abs(FTSignal/dots);
amplitude = amplitudedouble(1:dots/2+1);
amplitude(2:end-1) = 2*amplitude(2:end-1);
frequencies = (0:(dots/2))/(dots*period);

end
