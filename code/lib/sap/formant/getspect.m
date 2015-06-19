function [fftspect, lpcspect, freqscale, A, presig, poles] = getspect(frame, fs, lpcOrder)

% GETSPECT Get spectrum information and other LPC related information
%    The output results conclude FFT spectrum, LPC spectrum, and 
%    LPC parameters and predicted signal and poles.
%    Usage : [fftspect, lpcspect, freqscale, A, presig, poles] = getspect(frame, fs, lpcOrder).
%    
%    Cheng-Yuan Lin, 2003, January, 18.

fspace = 0: 1 : (fs/2);  %10 Hz spacing.

% Computer DFT, LPC, LPC spectrum
[freqscale, fftspect] = fftspectrum(frame, fs, 2, fspace);  %DFT
[freqscale, lpcspect, A] = lpcspectrum(frame, fs, lpcOrder, 2, fspace);  %LPC spectrum

% Compute predicted signal
A(1) = 0;
presig = filter(-A, 1, frame);

% Compute poles of LPC analysis.
A(1) = 1; 
poles = roots(A);

% Compute MSE as a percentage of signal MSE : this is reserved part.
%sigmse = sqrt(mean(seg.^2));
%errmse = sqrt(mean(err.^2));
%totalerror =  num2str(100*errmse/sigmse);
%errsig = frame - presig;
%errspectrum = performFFT(errsig, fspace, fs);
