function [freq, Pyy] = fftspectrum(frame, fs, fftmethod, fspace)
% FFTSPECTRUM Short term FFT to get magnitude spectrum information
%    Usage : [freq, Pyy] = fftspectrum(frame, fs, fftmethod, fspace)
%    frame means frame vector from some wavefile.
%    fs means sample rate of the wavefile.
%    fftmethod presents two different extraction rules. Its value is either 1 or 2.
%    fspace means frequency spacing.
%    freq represents frequency components
%    Pyy represents log magnitudes corresponding to frequency components.
%
%    Cheng-Yuan Lin 2003, January, 22.

if nargin<4, 
   fspace = 0: 20 : (fs/2);  %20 Hz spacing.
end;
if nargin<3, fftmethod = 1; end;

switch fftmethod,
case 1,
   % general version
   frameSize = length(frame);
   Y    = fft(frame);
   half = ceil(frameSize/2);
   Pyy  = 20*log10(abs(Y));
   Pyy  = Pyy(1:half)';  %fft spectrum is symmetric.
   %freq = fs/frameSize*(1:frameSize);
   %freq = freq(1:half);
   spectfactor = (fs/2)/length(Pyy);
   freq = [1: spectfactor : (fs/2)];
case 2,
   % Using interpolation version
   frameSize = length(frame);
   lf = length(fspace);
   f = fft(frame, max(128, power(2, round(log2(frameSize)))));
   f = 20*log10(1e-10 + abs(f(2:round(length(f)/2))));
   Pyy = interp1(linspace(0, fs/2, length(f)), f, fspace);
   freq = (fs/2)/length(Pyy)*(1:length(Pyy))';
end;