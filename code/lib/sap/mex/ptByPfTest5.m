% Compare wave2pitchByAcf.m and ptByPfMex.cpp
addMyPath;
clear all; close all;

fprintf('Compiling ptByPfMex.cpp...\n');
mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave ptByPfMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output ptByPfMex.dll

waveFile='../twinkle_twinkle_little_star.wav';
waveFile='../十個印第安人_不詳_0.wav';
[y, fs, nbits]=wavReadInt(waveFile);

% dataMode=1 ===> floating-point
% pfType=1 ===> acf
dataMode=1; pfType=1; ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.useVolThreshold=1; ptOpt.medianFilterOrder=0; ptOpt.cutLeadingTrailingZero=0;

tic
[pitch1, clarity1]=wave2pitchByAcf(y, ptOpt);
fprintf('Time for wave2pitchByAcf = %g seconds when pfType=%d and dataMode=%d\n', toc, pfType, dataMode);
tic
[pitch2, clarity2, pfMat]=ptByPfMex(y, fs, nbits, ptOpt, dataMode);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d and dataMode=%d\n', toc, pfType, dataMode);

subplot(4,1,1); plot((1:length(y))/fs, y); title(['wavFile=', waveFile]); axis tight
subplot(4,1,2); plot(1:length(pitch1), pitch1, '.-', 1:length(pitch2), pitch2/1024, '.-'); title('Pitch');
legend('Pitch 1 by wave2pitchByAcf', 'Pitch 2 by ptByPfmex');
subplot(4,1,3); plot(1:length(clarity1), clarity1, '.-', 1:length(clarity2), clarity2/1024, '.-'); title('Clarity');
legend('Clarity 1 by wave2pitchByAcf', 'Clarity 2 by ptByPfmex');
subplot(4,1,4); imagesc(flipud(pfMat));
