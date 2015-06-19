
fprintf('Compiling ptByDpOverPfMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility -I/users/jang/c/lib/wave ptByDpOverPfMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp /users/jang/c/lib/wave/waveRead4pda.cpp -output ptByDpOverPfMex
fprintf('Done compiling\n');

waveFile='但使龍城飛將在.wav';
[y, fs, nbits]=wavReadInt(waveFile);
pfType=1;	% 0 for AMDF, 1 for ACF
ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.useVolThreshold=0;
ptOpt.useClarityThreshold=0;
ptOpt.useSmooth=0;
[pitch, clarity, pfMat]=ptByDpOverPfMex(y, fs, nbits, ptOpt);

subplot(4,1,1); plot((1:length(y))/fs, y);
subplot(4,1,2); pitch(pitch==0)=nan; plot(pitch, '.-');
subplot(4,1,3); plot(clarity, '.-');
subplot(4,1,4); imagesc(pfMat); axis xy


figure;
%ptOpt=ptOptSet(8000, 8, pfType);
ptOpt.useEpd=1;
epdPrm=epdPrmSet(fs);
plotOpt=1;
[pitch, pitch2, pfMat]=ptByDpOverPf2(y, fs, nbits, ptOpt, epdPrm, plotOpt);