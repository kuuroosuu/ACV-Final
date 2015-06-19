mex epdByVolZcrMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

addMyPath
waveFile='../此恨綿綿無絕期+sil.wav';
%waveFile='C:\dataSet\資策會-派斌\pmp-epd\0001_16k.wav';
%waveFile='../csNthu_omniMic.wav';
[wave, fs, nbits]=wavReadInt(waveFile);
epdPrm.frameSize = round(fs/31.25);		% fs=8000 ===> frameSize=256;
epdPrm.overlap = round(epdPrm.frameSize/2);
epdPrm.volRatio1=0.1;		% lower threshold
epdPrm.volRatio2=0.2;		% upper threshold
epdPrm.zcrRatio=0.1;
epdPrm.zcrShiftGain=4;
plotOpt = 1;
ep1 = epdByVolZcr(wave, fs, nbits, epdPrm, plotOpt);
ep2 = epdByVolHodMex(wave, fs);
subplot(4,1,1);
line(ep2(1)*[1 1]/fs, [min(wave), max(wave)], 'color', 'r');
line(ep2(2)*[1 1]/fs, [min(wave), max(wave)], 'color', 'r');
if isequal(ep1, ep2)
	fprintf('Same endpoints!\n');
else
	fprintf('Different endpoints!\n');
	fprintf('Dev1: %g - %g = %g\n', ep1(1), ep2(1), ep1(1)-ep2(1));
	fprintf('Dev2: %g - %g = %g\n', ep1(2), ep2(2), ep1(2)-ep2(2));
end