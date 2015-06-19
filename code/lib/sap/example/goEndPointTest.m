clear mex
mex endpoint2mex.c

waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\coco1.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\soph5.wav';
[y, fs, nbits]=wavread(waveFile);
out2=endpoint2mex(y, fs, nbits)
out=endpoint2(y, fs, nbits)
difference=sum(abs(out-out2));
fprintf('difference=%d\n', difference);

return

waveDir='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機';
waveFiles=dir([waveDir, '/*.wav']);

for i=1:length(waveFiles)
	waveFile=[waveDir, '/', waveFiles(i).name];
	fprintf('%d/%d: %s ===> ', i, length(waveFiles), waveFile);
	[y, fs, nbits]=wavread(waveFile);
	out2=endpoint2mex(y, fs, nbits);
	out=endpoint2(y, fs, nbits);
	difference=sum(abs(out-out2));
	fprintf('difference=%d\n', difference);
	if difference~=0
		fprintf('Warning: difference is not zero!\n');
	end
end