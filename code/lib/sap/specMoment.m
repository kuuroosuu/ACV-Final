function [specMean, specStddev, specSkewness, specKurtosis]=specMoment(powerSpec, freq)
% specMoment: Compute the spectrum moments of power spectrum
%	Usage: [specMean, specStddev, specSkewness, specKurtosis]=specMoment(powerSpec, freq)
%		powerSpec: power spectrum (in linear scale)
%		freq: freq vector of the power spectrum
%		specMean, specStddev, specSkewness, specKurtosis: spectrum moments
%
%	For example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		[wave, fs, nbits]=wavReadInt(waveFile);
%		frameSize=256;
%		overlap=0;
%		frameMat=buffer2(wave, frameSize, overlap);
%		frameNum=size(frameMat, 2);
%		specMean=zeros(1, frameNum);
%		specStddev=zeros(1, frameNum);
%		specSkewness=zeros(1, frameNum);
%		specKurtosis=zeros(1, frameNum);
%		powerDbMat=[];
%		for i=1:frameNum
%			[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
%			powerDbMat=[powerDbMat, powerDb];
%			powerSpec=magSpec.^2;
%			[specMean(i), specStddev(i), specSkewness(i), specKurtosis(i)]=specMoment(powerSpec, freq);
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
%		moments=[specMean; specStddev; specSkewness; specKurtosis]';
%		plot(frameTime, moments);
%		set(gca, 'xlim', [min(frameTime), max(frameTime)]);
%		xlabel('Time (sec)'); ylabel('Moments');
%		legend('Mean', 'Stddev', 'Skewness', 'Kurtosis');

%	Roger Jang, 20070506, 20080401

if nargin<1, selfdemo; return; end

nps=powerSpec/sum(powerSpec);		% Normalized power spectrum
specMean=dot(freq, nps);
specStddev=sqrt(dot((freq-specMean).^2, nps));
specSkewness=dot(((freq-specMean)/specStddev).^3, nps);
specKurtosis=-3+dot(((freq-specMean)/specStddev).^4, nps);

% ====== Self demo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
[wave, fs, nbits]=wavReadInt(waveFile);
frameSize=256;
overlap=0;
frameMat=buffer2(wave, frameSize, overlap);
frameNum=size(frameMat, 2);
specMean=zeros(1, frameNum);
specStddev=zeros(1, frameNum);
specSkewness=zeros(1, frameNum);
specKurtosis=zeros(1, frameNum);
powerDbMat=[];
for i=1:frameNum
	[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
	powerDbMat=[powerDbMat, powerDb];
	powerSpec=magSpec.^2;
	[specMean(i), specStddev(i), specSkewness(i), specKurtosis(i)]=specMoment(powerSpec, freq);
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
moments=[specMean; specStddev; specSkewness; specKurtosis]';
plot(frameTime, moments);
set(gca, 'xlim', [min(frameTime), max(frameTime)]);
xlabel('Time (sec)'); ylabel('Moments');
legend('Mean', 'Stddev', 'Skewness', 'Kurtosis');