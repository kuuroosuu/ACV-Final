waveFile='zh_ch_sh_z_c_s.wav';
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
freqBandEnergyRatio=sum(freqBandEnergy(2:3, :))./freqBandEnergy(1, :);
plot(frameTime, freqBandEnergyRatio');
set(gca, 'xlim', [min(frameTime), max(frameTime)]);
xlabel('Time (sec)'); ylabel('Frequency band energy ratio');
%legend('0~fs/6', 'fs/6~fs/3', 'fs/3~fs/2');