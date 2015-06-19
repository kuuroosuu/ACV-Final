function rolloff=specRolloff(powerSpec, freq, ratio)
% specRolloff: Compute the rolloff frequency of a given power spectrum and a given ratio
%	Usage: roffoff=specRolloff(powerSpec, freq, ratio)
%		powerSpec: spectrum (in linear scale)
%		freq: freq vector of the power spectrum
%		ratio: a given ratio
%		roffoff: roffoff frequency
%
%	For example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		[wave, fs, nbits]=wavReadInt(waveFile);
%		frameSize=256;
%		overlap=0;
%		frameMat=buffer2(wave, frameSize, overlap);
%		frameNum=size(frameMat, 2);
%		rolloff=zeros(1, frameNum);
%		powerDbMat=[];
%		for i=1:frameNum
%			[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
%			powerDbMat=[powerDbMat, powerDb];
%			powerSpec=magSpec.^2;
%			rolloff(i)=specRolloff(powerSpec, freq);
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
%		plot(frameTime, rolloff, '.-');
%		set(gca, 'xlim', [min(frameTime), max(frameTime)]);
%		xlabel('Time (sec)'); ylabel('Rolloff frequency (Hz)');

%	Roger Jang, 20070506, 20080401

if nargin<1, selfdemo; return; end
if nargin<3, ratio=0.75; end

total=sum(powerSpec);
threshold=total*ratio;
partial=0;
for i=1:length(powerSpec)
	partial=partial+powerSpec(i);
	if (partial>=threshold)
		rolloff=freq(i);
		return;
	end
end
rolloff=freq(end);

% ====== Self demo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
[wave, fs, nbits]=wavReadInt(waveFile);
frameSize=256;
overlap=0;
frameMat=buffer2(wave, frameSize, overlap);
frameNum=size(frameMat, 2);
rolloff=zeros(1, frameNum);
powerDbMat=[];
for i=1:frameNum
	[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
	powerDbMat=[powerDbMat, powerDb];
	powerSpec=magSpec.^2;
	rolloff(i)=specRolloff(powerSpec, freq);
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
plot(frameTime, rolloff, '.-');
set(gca, 'xlim', [min(frameTime), max(frameTime)]);
xlabel('Time (sec)'); ylabel('Rolloff frequency (Hz)');