function freqBandEnergy=specFreqBandEnergy(powerSpec, freq, freqBand)
% specMoment: Compute the frequency band energy of power spectrum
%	Usage: freqBandEnergy=specFreqBandEnergy(powerSpec, freq, freqBand)
%		powerSpec: power spectrum (in linear scale)
%		freq: freq vector of the power spectrum
%		freqBand: freq bands in the format of [L1, U1, L2, U2, L3, U3, ...]
%
%	For example:
%		waveFile='zh_z_ch_c_sh_s.wav';
%		[wave, fs, nbits]=wavReadInt(waveFile);
%		frameSize=256;
%		overlap=0;
%		frameMat=buffer2(wave, frameSize, overlap);
%		frameNum=size(frameMat, 2);
%		maxFreq=fs/2;
%		freqBand=[0, maxFreq/3, maxFreq/3+eps, maxFreq*2/3, maxFreq*2/3+eps, maxFreq];
%		freqBandEnerty=zeros(length(freqBand)/2, frameNum);
%		powerDbMat=[];
%		for i=1:frameNum
%			[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
%			powerDbMat=[powerDbMat, powerDb];
%			powerSpec=magSpec.^2;
%			freqBandEnergy(:,i)=specFreqBandEnergy(powerSpec, freq, freqBand);
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
%		plot(frameTime, freqBandEnergy');
%		set(gca, 'xlim', [min(frameTime), max(frameTime)]);
%		xlabel('Time (sec)'); ylabel('Frequency band energy');
%		legend('0~fs/6', 'fs/6~fs/3', 'fs/3~fs/2');

%	Roger Jang, 20080401

if nargin<1, selfdemo; return; end

freqBandNum=length(freqBand)/2;
freqBandEnergy=zeros(freqBandNum, 1);
for i=1:freqBandNum
	lower=freqBand(2*i-1);
	upper=freqBand(2*i);
	freqBandEnergy(i, 1)=sum(powerSpec(lower<=freq & freq<upper));
end
	
% ====== Self demo
function selfdemo
waveFile='zh_z_ch_c_sh_s.wav';
[wave, fs, nbits]=wavReadInt(waveFile);
frameSize=256;
overlap=0;
frameMat=buffer2(wave, frameSize, overlap);
frameNum=size(frameMat, 2);
maxFreq=fs/2;
freqBand=[0, maxFreq/3, maxFreq/3+eps, maxFreq*2/3, maxFreq*2/3+eps, maxFreq];
freqBandEnerty=zeros(length(freqBand)/2, frameNum);
powerDbMat=[];
for i=1:frameNum
	[magSpec, phaseSpec, freq, powerDb]=fftOneSide(frameMat(:,i), fs);
	powerDbMat=[powerDbMat, powerDb];
	powerSpec=magSpec.^2;
	freqBandEnergy(:,i)=specFreqBandEnergy(powerSpec, freq, freqBand);
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
plot(frameTime, freqBandEnergy');
set(gca, 'xlim', [min(frameTime), max(frameTime)]);
xlabel('Time (sec)'); ylabel('Frequency band energy');
legend('0~fs/6', 'fs/6~fs/3', 'fs/3~fs/2');