function [mag, phase, freq, powerDb]=fftTwoSide(signal, fs, plotOpt)
% fftTwoSide: Two-sided FFT for real/complex signals
%	Usage: [mag, phase, freq, powerDb]=fftTwoSide(signal, fs)

%	Roger Jang, 20060411

if nargin<1, selfdemo; return; end
if nargin<2, fs=1; end
if nargin<3, plotOpt=0; end

N = length(signal);		% 翴计
freqStep = fs/N;		% 繵办繵瞯秆猂
time = (0:N-1)/fs;		% 办丁ㄨ
z = fft(signal);		% Spectrum
freq = freqStep*(-N/2:N/2-1);	% 繵办繵瞯ㄨ
z = fftshift(z);		% 盢繵瞯禸箂翴竚い
mag=abs(z);			% Magnitude
phase=unwrap(angle(z));		% Phase
powerDb=20*log(mag+eps);	% Power in db

if plotOpt
	% ====== Plot time-domain signals
	subplot(3,1,1);
	plot(time, signal, '.-');
	title(sprintf('Input signals (fs=%d)', fs));
	xlabel('Time (seconds)'); ylabel('Amplitude'); axis tight
	% ====== Plot spectral power
	subplot(3,1,2);
	plot(freq, powerDb, '.-'); grid on
	title('Power spectrum');
	xlabel('Frequency (Hz)'); ylabel('Power (db)'); axis tight
	% ====== Plot phase
	subplot(3,1,3);
	plot(freq, phase, '.-'); grid on
	title('Phase');
	xlabel('Frequency (Hz)'); ylabel('Phase (Radian)'); axis tight
end

% ====== Self demo
function selfdemo
[y, fs]=wavread('welcome.wav');
frameSize=512;
startIndex=2047;
signal=y(startIndex:startIndex+frameSize+1);
signal=signal.*hamming(length(signal));
feval(mfilename, signal, fs, 1);