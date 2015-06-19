waveFile='speechWith50HzNoise02.wav';
[y, fs, nbits]=wavReadInt(waveFile);
frameSize=320;
overlap=0;
frameMat=buffer2(y, frameSize, overlap);
frameNum=size(frameMat, 2);
frame=frameMat(:, round(frameNum/2+3));

subplot(2,1,1);
time=(1:length(y))/fs;
plot(time, y);
set(gca, 'xlim', [min(time), max(time)]);
xlabel('Time (sec)'); title(sprintf('Waveform of "%s"', strrep(waveFile, '_', '\_')));

subplot(2,1,2);
plot(frame);

[mag, phase, freq, powerSpec]=fftOneSide(y, fs, 1);
