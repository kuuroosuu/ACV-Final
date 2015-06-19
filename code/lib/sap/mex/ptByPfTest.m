addMyPath;
%fprintf('Compiling ptByPfMex.cpp...\n');
%mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave ptByPfMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output ptByPfMex.dll

waveFile='..\主人下馬客在船.wav';
waveFile='..\楊家有女初長成.wav';
waveFile='..\但使龍城飛將在.wav';
waveFile='..\twinkle_twinkle_little_star.wav';

[y, fs, nbits]=wavReadInt(waveFile);
fprintf('waveFile=%s\n', waveFile);

tic
pfType=0; ptOpt=ptOptSet(fs, nbits, pfType); ptOpt.cutLeadingTrailingZero=0;
[pitch1, clarity1]=ptByPfMex(y, fs, nbits, ptOpt, 1);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d (AMDF)\n', toc, pfType);

tic
pfType=1; ptOpt=ptOptSet(fs, nbits, pfType); ptOpt.cutLeadingTrailingZero=0;
[pitch2, clarity2]=ptByPfMex(y, fs, nbits, ptOpt, 1);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d (ACF)\n', toc, pfType);

tic
pfType=1; ptOpt=ptOptSet(fs, nbits, pfType); ptOpt.cutLeadingTrailingZero=0;
ptOpt.dpUseLocalOptim=1;
ptOpt.pfWeight=1;
ptOpt.indexDiffWeight=100;
pitch3=1024*ptByDpOverPfMex(y, fs, nbits, ptOpt);
fprintf('Time for ptByDpOverPfMex = %g seconds\n', toc);

subplot(3,1,1);
plot((1:length(y))/fs, y);
subplot(3,1,2);
pitch1(pitch1==0)=nan; pitch2(pitch2==0)=nan; pitch3(pitch3==0)=nan;
plot(1:length(pitch1), pitch1, '.-', 1:length(pitch2), pitch2, '.-', 1:length(pitch3), pitch3, '.-');
ylabel('Pitch'); legend('pfType=1 (AMDF)', 'pfType=2 (ACF)', 'DP');
subplot(3,1,3);
plot(1:length(clarity1), clarity1, '.-', 1:length(clarity2), clarity2, '.-');
ylabel('Clarity'); legend('pfType=1 (AMDF)', 'pfType=2 (ACF)');