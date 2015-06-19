waveFile='\dataSet\dean\201005-mirRecording\Record_Wave/陳亮宇_1/中華民國國歌_不詳_0_s.wav';
[y, fs, nbits]=wavread(waveFile);
frameSize=512;
overlap=0;
frameMat=buffer2(y, frameSize, overlap);
frame=frameMat(:,54);
[magSpec, phaseSpec, freq, ps]=fftOneSide(frame, fs);
n=length(ps);
order=10;
ave=polyval(polyfit((1:n)', ps, order), (1:n)');
subplot(4,1,1);
plot(1:n, ps, 1:n, ave);
subplot(4,1,2);
plot(1:n, ps-ave);

subplot(4,1,3);
amdf=frame2amdfMex(ps-ave);
plot(amdf);

subplot(4,1,4);
acf=frame2nsdfMex(ps-ave);
plot(acf);

figure;
nsdf=frame;
for i=1:4
	subplot(4,1,i);
	plot(nsdf);
	nsdf=frame2acfMex(nsdf, length(nsdf), 2);
end
