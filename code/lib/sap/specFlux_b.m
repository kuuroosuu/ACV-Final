function flux=specFlux(powerSpec, freq)
% specFlux: Compute the flux of power spectrum
%	Usage: flux=specFlux(powerSpec, freq, freqRange)
%		powerSpec: power spectrum (in linear scale)
%		freq: freq vector of the power spectrum
%		flux: output flux
%
%	For example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		[wave, fs, nbits]=wavReadInt(waveFile);
%		frameSize=256;
%		overlap=0;
%		frameMat=buffer2(wave, frameSize, overlap);
%		frameNum=size(frameMat, 2);
%		flux=zeros(1, frameNum);
%		powerDbMat=[];
%		for i=1:frameNum
%			[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
%			powerDbMat=[powerDbMat, powerDb];
%			powerSpec=magSpec.^2;
%			flux(i)=specFlux(powerSpec, freq);
%		end
%		subplot(3,1,1);
%		time=(1:length(wave))/fs;
%		plot(time, wave);
%		set(gca, 'xlim', [min(time), max(time)]);
%		xlabel('Time (sec)'); title(sprintf('Waveform of "%s"', strrep(waveFile, '_', '\_')));
%		subplot(3,1,2);
%		imagesc(powerDbMat); axis xy
%		title('Spectrogram');
%		subplot(3,1,3);
%		frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
%		plot(frameTime, flux, '.-');
%		set(gca, 'xlim', [min(frameTime), max(frameTime)]);
%		xlabel('Time (sec)'); ylabel('Flux');

%	Roger Jang, 20070506, 20080401

if nargin<1, selfdemo; return; end

flux=sum(abs(powerSpec(1:end-1)-powerSpec(2:end)));

% ====== Self demo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
[wave, fs, nbits]=wavReadInt(waveFile);
frameSize=256;
overlap=0;
frameMat=buffer2(wave, frameSize, overlap);
frameNum=size(frameMat, 2);
flux=zeros(1, frameNum);
powerDbMat=[];
for i=1:frameNum
	[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
	powerDbMat=[powerDbMat, powerDb];
	powerSpec=magSpec.^2;
	flux(i)=specFlux(powerSpec, freq);
end
subplot(3,1,1);
time=(1:length(wave))/fs;
plot(time, wave);
set(gca, 'xlim', [min(time), max(time)]);
xlabel('Time (sec)'); title(sprintf('Waveform of "%s"', strrep(waveFile, '_', '\_')));
subplot(3,1,2);
imagesc(powerDbMat); axis xy
title('Spectrogram');
subplot(3,1,3);
frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
plot(frameTime, flux, '.-');
set(gca, 'xlim', [min(frameTime), max(frameTime)]);
xlabel('Time (sec)'); ylabel('Flux');