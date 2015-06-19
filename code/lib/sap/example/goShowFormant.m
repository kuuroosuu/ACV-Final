waveFile='但使龍城飛將在.wav';
formantPrm.formantNum=3;		% No. of formants
formantPrm.frameSize=20;		% 20 ms
formantPrm.frameStep=10;		% 10 ms
formantPrm.lpcOrder=12;		% order of LPC
formant=wave2formant(waveFile, formantPrm, 1);

frameMat=buffer2(y, frameSize, overlap);
frameNum=size(frameMat, 2);
plot(formant);

return

[y, fs, nbits]=wavReadInt(waveFile);
frameSize=formantPrm.frameSize*fs/1000;
frameStep=formantPrm.frameStep*fs/1000;
overlap=frameSize-frameStep;
frameMat=buffer2(y, frameSize, overlap);
spectrogram(y, frameSize, frameStep, 512, fs);