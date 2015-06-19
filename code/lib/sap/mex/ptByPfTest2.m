% Check the difference between different dataMode of int and double
addMyPath;
clear all; close all;

%fprintf('Compiling ptByPfMex.cpp...\n');
%mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave ptByPfMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output ptByPfMex.dll

waveFile='..\主人下馬客在船.wav';
waveFile='..\楊家有女初長成.wav';
waveFile='..\但使龍城飛將在.wav';
waveFile='..\twinkle_twinkle_little_star.wav';

[y, fs, nbits]=wavReadInt(waveFile);
fprintf('waveFile=%s\n', waveFile);

pfType=1; ptOpt=ptOptSet(fs, nbits, pfType); ptOpt.medianFilterOrder=0; ptOpt.cutLeadingTrailingZero=0;

tic
dataMode=0;	% int data type
[pitch1, clarity1]=ptByPfMex(y, fs, nbits, ptOpt, dataMode);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d and dataMode=%d\n', toc, pfType, dataMode);

tic
dataMode=1;	% float data type
[pitch2, clarity2]=ptByPfMex(y, fs, nbits, ptOpt, dataMode);
fprintf('Time for ptByPfMex = %g seconds when pfType=%d and dataMode=%d\n', toc, pfType, dataMode);

subplot(2,2,1); plot([pitch1', pitch2'], '.-'); title('Pitch of different data modes'); legend('pitch1', 'pitch2');
subplot(2,2,2); plot([clarity1', clarity2'], '.-'); title('Clarity of different data modes'); legend('clarity1', 'clarity2');
subplot(2,2,3); plot(abs(pitch1-pitch2), '.-'); title('Diff. in pitch of different data modes');
subplot(2,2,4); plot(abs(clarity1-clarity2), '.-'); title('Diff. in clarity of different data modes');
