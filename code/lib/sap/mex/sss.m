addMyPath;

fprintf('Compiling ptByPfMex.cpp...\n');
mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave ptByPfMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output ptByPfMex.dll

waveFile='25 Minutes_M.L.T.R_0.wav';

[y, fs, nbits]=wavReadInt(waveFile);
fprintf('waveFile=%s\n', waveFile);

tic
pfType=2; ptOpt=ptOptSet(fs, nbits, pfType);
[pitch2, clarity2]=ptByPfMex(y, fs, nbits, ptOpt);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d (ACF)\n', toc, pfType);
