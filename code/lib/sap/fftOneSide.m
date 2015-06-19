function [magSpec, phaseSpec, freq, powerSpecInDb]=fftOneSide(frameMat, fs, plotOpt)
% fftOneSide: One-sided FFT for real signals
%	Usage: [magSpec, phaseSpec, freq, powerSpecInDb]=fftOneSide(frameMat, fs, plotOpt)
%
%	For example:
%		[y, fs]=wavread('welcome.wav');
%		frameSize=512;
%		startIndex=2047;
%		frame=y(startIndex:startIndex+frameSize+1);
%		frame=frame.*hamming(length(frame));
%		plotOpt=1;
%		[magSpec, phaseSpec, freq, powerSpecInDb]=fftOneSide(frame, fs, plotOpt);

%	Roger Jang, 20060411, 20070506

if nargin<1, selfdemo; return; end
if nargin<2, fs=1; end
if nargin<3, plotOpt=0; end
[dim, frameCount]=size(frameMat);

if isstruct(fs)		% fs is in fact prm
	prm=fs;
	if ~isfield(prm, 'plotOpt'), prm.plotOpt=0; end
	if ~isfield(prm, 'nfft'), prm.nfft=dim; end
elseif isempty(fs)	% Set default prm
	prm.fs=1;
	prm.plotOpt=0;
	prm.nfft=dim;
else				% Create prm for each input arguments
	prm.fs=fs;
	prm.plotOpt=plotOpt;
	prm.nfft=dim;
end
fs=prm.fs;
plotOpt=prm.plotOpt;
nfft=prm.nfft;

[frameSize, frameNum]=size(frameMat);
if frameSize==1, frameMat=frameMat(:); end	% Single frame as a row vector

freqStep = fs/frameSize;		% Frequency resolution
time = (0:frameSize-1)/fs;		% Time vector

if nfft>dim
	frameMat=[frameMat; zeros(nfft-dim, frameNum)];
end
z = fft(frameMat);			% Spectrum
freq = freqStep*(0:frameSize/2)';	% Frequency vector
z = z(1:length(freq), :);		% One side
z(2:end-1, :)=2*z(2:end-1, :);		% Assuming N is even, symmetric data is multiplied by 2
magSpec=abs(z);				% Magnitude spectrum
phaseSpec=unwrap(angle(z));		% Phase spectrum
%powerSpecInDb=20*log(magSpec+realmin);	% Power in db
powerSpecInDb=20*log(magSpec);		% Power in db

if plotOpt & frameNum==1
	% ====== Plot time-domain signals
	subplot(3,1,1);
	plot(time, frameMat, '.-');
	title(sprintf('Input signals (fs=%d)', fs));
	xlabel('Time (seconds)'); ylabel('Amplitude'); axis tight
	% ====== Plot spectral power
	subplot(3,1,2);
	plot(freq, powerSpecInDb, '.-'); grid on
	title('Power spectrum');
	xlabel('Frequency (Hz)'); ylabel('Power (db)'); axis tight
	% ====== Plot phase
	subplot(3,1,3);
	plot(freq, phaseSpec, '.-'); grid on
	title('Phase');
	xlabel('Frequency (Hz)'); ylabel('Phase (Radian)'); axis tight
end

% ====== Self demo
function selfdemo
[y, fs]=wavread('welcome.wav');
frameSize=512;
startIndex=2047;
frame=y(startIndex:startIndex+frameSize+1);
frame=frame.*hamming(length(frame));
%frame=[frame; zeros(frameSize, 1)];
[magSpec, phaseSpec, freq, powerSpecInDb]=feval(mfilename, frame, fs, 1);