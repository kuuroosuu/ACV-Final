addMyPath;

clear all;
fprintf('Compiling wave2pfMatMex.cpp...\n');
%mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave wave2pfMatMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp -output wave2pfMatMex.dll
mex -Id:/users/jang/c/lib -Id:/users/jang/c/lib/audio -Id:/users/jang/c/lib/utility -Id:/users/jang/c/lib/wave wave2pfMatMex.cpp d:/users/jang/c/lib/audio/audio.cpp d:/users/jang/c/lib/utility/utility.cpp d:/users/jang/c/lib/wave/waveRead4pda.cpp
fprintf('Done with compiling\n');
dos('copy /y wave2pfMatMex.mex* ..\private');

waveFile='yankeeDoodle_8k8b.wav';
waveFile='ivegotrhythm_16k16b.wav';
waveFile='..\山月照彈琴.wav';
waveFile='..\但使龍城飛將在.wav';
fprintf('waveFile=%s\n', waveFile);

[y, fs, nbits]=wavRead(waveFile);
ptOpt=ptOptSet(fs, nbits);
pfType=0;	% 0 for AMDF, 1 for ACF
if pfType==0, pfMethod=2; end
if pfType==1, pfMethod=1; end

tic, pfMat1 = wave2pfMatM  (y, ptOpt.frameSize, ptOpt.overlap, ptOpt.pfLen, pfType, pfMethod); time1=toc;
tic, pfMat2 = wave2pfMatMex(y, ptOpt.frameSize, ptOpt.overlap, ptOpt.pfLen, pfType, pfMethod); time2=toc;

subplot(2,1,1);pcolor(pfMat1);  shading interp
subplot(2,1,2);pcolor(pfMat2); shading interp

%pfMat=pfMat(:,1);
%pfMat2=pfMat2(:,1);
figure;
subplot(3,1,1); plot(1:ptOpt.pfLen, pfMat1);
subplot(3,1,2); plot(1:ptOpt.pfLen, pfMat2);
subplot(3,1,3); plot(1:ptOpt.pfLen, pfMat1-pfMat2);

fprintf('difference in pfMat = %d\n', sum(sum(abs(pfMat1-pfMat2))));
fprintf('Time for m = %f\n', time1);
fprintf('Time for mex = %f\n', time2);
fprintf('Time ratio = %f\n', time1/time2);
