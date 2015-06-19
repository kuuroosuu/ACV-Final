%mex endPointDetectMex.cpp audio.cpp d:\users\jang\c\lib\lib.cpp
mex endPointDetectMex.cpp

addMyPath
waveFile='¦¹«ëºøºøµLµ´´Á+sil.wav';
[y, fs, nbits]=wavreadInt(waveFile);
[ep1, vol1, zcr1] = endPointDetect(y, fs, nbits);
%[ep2, vol2, zcr2] = endPointDetectMex(y, fs, nbits);
ep2 = endPointDetectMex(y, fs, nbits);
if isequal(ep1, ep2)
	fprintf('Same endpoints!\n');
end
return

frameSize=320;
overlap=240;
framedY=buffer2(y, frameSize, overlap);
volume=frame2volume(framedY);
frameNum=size(framedY, 2);
time=(1:length(y))/fs;
frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;

subplot(3,1,1);
plot(time, y); axis tight
line(ep1(1)/fs*[1 1], [-1 1]*2^nbits/2, 'color', 'r');
line(ep1(2)/fs*[1 1], [-1 1]*2^nbits/2, 'color', 'r');
line(ep2(1)/fs*[1 1], [-1 1]*2^nbits/2, 'color', 'm');
line(ep2(2)/fs*[1 1], [-1 1]*2^nbits/2, 'color', 'm');
subplot(3,1,2);
plot(frameTime, vol1, frameTime, vol2); axis tight
fprintf('Volume difference = %g\n', sum(abs(vol1-vol2)));
subplot(3,1,3);
plot(frameTime, zcr1, frameTime, zcr2); axis tight
fprintf('Zcr difference = %g\n', sum(abs(zcr1-zcr2)));