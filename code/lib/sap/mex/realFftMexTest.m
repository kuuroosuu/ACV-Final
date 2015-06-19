waveFile='yankeeDoodle_8k8b.wav';
[y, fs, nbits]=wavread(waveFile);
frameSize=256;
overlap=0;
useHammingWin=1;

win=hamming(frameSize);	% Make it compatible with specgram
tic
specMat1 = abs(spectrogram(y, win, 0, frameSize, fs));
fprintf('Time for spectrogram = %.2g sec\n', toc)
frameMat=buffer2(y, frameSize, overlap);
tic
specMat2=realFftMex(frameMat, useHammingWin);
fprintf('Time for realFftMex = %.2g sec\n', toc)
specMat2=specMat2(1:frameSize/2+1, :);

subplot(211); imagesc(specMat1); axis xy
subplot(212); imagesc(specMat2); axis xy

fprintf('isequal(specMat1, specMat2)=%d\n', isequal(specMat1, specMat2));
fprintf('Max deviation = %f\n', max(max(abs(specMat1-specMat2))));

