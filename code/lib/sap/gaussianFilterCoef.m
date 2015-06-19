function b=gaussianFilterCoef(fs, fco, plotOpt)
% gaussianFilterCoef: Return coefficients of Gaussian lowpass filter.
%
%	Usage:
%		b=gaussianFilterCoef(fs, fco, plotOpt)
%
%	Description:
%		b=gaussianFilterCoef(fs, fco) returns the coef. of Gaussian lowpass filter
%			fs=sampling rate
%			fco=cutoff (-3dB) freq, in Hz. 
%		Coeffs for FIR filter of length L (L always odd) are computed.
%		This symmetric FIR filter of length L=2N+1 has delay N/fs seconds.
%		(This function is modified from Wililiam Rose's code of gaussfiltcoef.m.)
%
%	Example:
%		% === Plot the freq response
%		fs=8000;
%		fco=300;
%		b=gaussianFilterCoef(fs, fco, 1);
%		% === Plot the signals after the filter
%		time=(1:100)/fs;
%		x=randn(length(time), 1);
%		x=sin(2*pi*200*time);
%		y=filter(b, 1, x);
%		figure; plot(time, x, time, y);
%		xlabel('Time (sec)'); ylabel('Amplitude');
%		legend('Original', 'After filtering');

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

b=0;
a=3.011*fco;
N=ceil(0.398*fs/fco);   %filter half-width, excluding midpoint
%Width N corresponds to at least +-3 sigma which captures at least 99.75%
%of area under Normal density function. sigma=1/(a*sqrt(2pi)).
L=2*N+1;                %full length of FIR filter
for k=-N:N
    b(k+N+1)=3.011*(fco/fs)*exp(-pi*(a*k/fs)^2);
end;
%b(k) coeffs computed above will add to almost exactly unity, but not 
%quite exact due to finite sampling and truncation at +- 3 sigma.  
%Next line adjusts to make coeffs b(k) sum to exactly unity.
b=b/sum(b);

if plotOpt
	subplot(2,1,1); plot(b); title('Magnitude of coefficients "b"');
	subplot(2,1,2);
	[h, w]=freqz(b, 1, 256, fs);
	plot(w, 10*log(abs(h))); xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
	title('Freq. response of the filter');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
