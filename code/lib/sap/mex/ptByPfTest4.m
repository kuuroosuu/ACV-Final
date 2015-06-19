addMyPath;
clear all; close all;

%fprintf('Compiling ptByPfMex.cpp...\n');
%mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave ptByPfMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output ptByPfMex.dll

waveFile='You Are My Sunshine_คฃธิ_0##.wav';
waveFile='yankeeDoodle_8k8b.wav';

[y, fs, nbits]=wavReadInt(waveFile);
%y=y/256;
fprintf('waveFile=%s\n', waveFile);

tic
% dataMode=1 ===> floating-point
% pfType=1 ===> acf
dataMode=1; pfType=1; ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.useVolThreshold=0; ptOpt.medianFilterOrder=0; ptOpt.cutLeadingTrailingZero=0;
[pitch1, clarity1, pfMat]=ptByPfMex(y, fs, nbits, ptOpt, dataMode);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d and dataMode=%d\n', toc, pfType, dataMode);

subplot(4,1,1); plot((1:length(y))/fs, y); axis tight
subplot(4,1,2); plot(pitch1, '.-'); title('Pitch');
subplot(4,1,3); plot(clarity1, '.-'); title('Clarity');
subplot(4,1,4); imagesc(flipud(pfMat));
