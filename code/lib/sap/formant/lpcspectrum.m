function [freq, Pyy, A] = lpcspectrum(frame, fs, lpcOrder, lpcmethod, fspace)
% LPCSPECTRUM Short term LPC to get magnitude spectrum information
%    Usage : [freq, Pyy, A] = lpcspectrum(frame, fs, lpcOrder, lpcmethod, fspace)
%    frame means frame vector from some wavefile.
%    fs means sample rate of the wavefile.
%    lpcOrder means order of LPC coefficients.
%    lpcmethod presents two different extraction rules. Its value is either 1 or 2.
%    fspace means frequency spacing.
%    freq represents frequency components
%    Pyy represents log magnitudes corresponding to frequency components.
%    A means LPC coefficients.
%
%    Cheng-Yuan Lin 2003, January, 22.

if nargin<5, 
   fspace = 0: 10 : (fs/2);  %10 Hz spacing.
end;
if nargin<4, lpcmethod = 2; end;

switch lpcmethod,
case 1,
   % general version
   frameSize = length(frame);
   [A, K, E] = ldlpc(frame, lpcOrder);
   sequence = [A zeros(1, frameSize-lpcOrder-1)];
   half = ceil(frameSize/2);
   Pyy  = 1./abs(fft(sequence));
   Pyy  = Pyy(1:half);  %fft spectrum is symmetric.
   Pyy  = 20*log10(Pyy);
   spectfactor = (fs/2)/length(Pyy);
   freq = [1: spectfactor : (fs/2)];
case 2,
   % Using MATLAB function : freqz.
   A    = real(lpc(frame, lpcOrder));
   Pyy  = real(20*log10(freqz(1, A, fspace, fs)));  % LPC spectrum
   freq = (fs/2)/length(Pyy)*(1:length(Pyy))';
end;