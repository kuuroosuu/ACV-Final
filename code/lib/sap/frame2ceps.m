function ceps=frame2ceps(frame, fs, zeroPaddedFactor, plotOpt)
% frame2hps: Conversion from a frame to harmonic product spectrum, for pitch tracking
%
%	Usage:
%		out=hps(frame, zeroPaddedFactor, plotOpt)
%
%	Example:
%		waveFile='soo.wav';
%		[y, fs, nbits] = wavread(waveFile);
%		startIndex=15000;
%		frameSize=round(32*fs/1000);
%		zeroPaddedFactor=15	% zero padding
%		frame=y(startIndex:startIndex+frameSize-1);
%		ceps=frame2ceps(frame, fs, zeroPaddedFactor, 1);

%	Roger Jang, 20080520

if nargin<1; selfdemo; return; end
if nargin<3, zeroPaddedFactor=15; end
if nargin<4; plotOpt=0; end

frame2=frame.*hamming(length(frame));	% Hamming windowing
frameSize=length(frame);
frame2=[frame2; zeros(zeroPaddedFactor*frameSize, 1)];		% zero padding
fftResult=log(abs(fft(frame2)));
%[mag, phase, freq, fftResult]=fftOneSide(frame2, fs);
ceps=ifft(fftResult);

if plotOpt
	plotNum=3;
	subplot(plotNum,1,1);
	plot(1:length(frame), frame, '.-'); axis tight;
	title('Frame');
	subplot(plotNum,1,2);
	plot(1:length(fftResult), fftResult);
	axis tight;
	title('Power spectrum');	
	subplot(plotNum,1,3);
	plot(1:length(ceps), ceps);
	axis tight;
	axisLimit=axis;
	title('Cepstrum');
	% ====== Compute pitch
	maxFreq=1000;
	minFreq=40;
	n1=round(fs/maxFreq);	% pdf(1:n1) will not be used
	n2=round(fs/minFreq);	% pdf(n2:end) will not be used
	ceps2=ceps; ceps2(1:n1)=-inf; ceps2(n2:end)=-inf;
	[maxValue, maxIndex]=max(ceps2);
	pitch=freq2pitch(fs/(maxIndex-1));
	line(n1*[1 1], axisLimit(3:4), 'color', 'r');
	line(n2*[1 1], axisLimit(3:4), 'color', 'r');
	line(maxIndex, ceps(maxIndex), 'marker', 'o', 'color', 'r');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
