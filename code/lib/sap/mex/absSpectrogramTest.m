fprintf('Compiling absSpectrogramMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility absSpectrogramMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp

fprintf('Testing absSpectrogramMex.cpp...\n');
waveFile='yankeeDoodle_8k8b.wav';
[y, fs, nbits]=wavread(waveFile);
frameSize=256;
overlap=0;
paddedSize=512;
useHammingWin=1;

win=hamming(frameSize);	% Make it compatible with specgram
tic
specMat1=abs(spectrogram(y, win, overlap, paddedSize, fs));
fprintf('Time for spectrogram = %.2g sec\n', toc)
tic
specMat2=absSpectrogramMex(y, frameSize, overlap, paddedSize, useHammingWin);
fprintf('Time for absSpectrogramMex = %.2g sec\n', toc)

subplot(211); imagesc(specMat1); axis xy
subplot(212); imagesc(specMat2); axis xy

fprintf('isequal(specMat1, specMat2)=%d\n', isequal(specMat1, specMat2));
maxDiff=max(max(abs(specMat1-specMat2)));
fprintf('Max deviation = %g\n', maxDiff);
