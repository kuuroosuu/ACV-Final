close all
waveFile='18.wav';
[y, fs, nbits]=wavread(waveFile);
y=y(:,1);
frameSize=256;
overlap=0;
framedY=buffer(y, frameSize, overlap);
volume=sum(abs(framedY));
plot(volume);

[segment, threshold]=splitVector(volume, 10*fs/frameSize, 0.5*fs/frameSize, 1);

for i=1:length(segment)
	sampleIndex1=(segment(i).beginIndex-1)*frameSize+1;
	sampleIndex2=(segment(i).endIndex)*frameSize;
	fprintf('Hit return to hear segment %d ...\n', i); pause
	sound(y(sampleIndex1:sampleIndex2), fs);
end