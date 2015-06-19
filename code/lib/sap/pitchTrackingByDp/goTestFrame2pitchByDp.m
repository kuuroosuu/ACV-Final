% 測試 weight 對於音高曲線的影響

waveFile='mountain.wav';
[y, fs, nbits]=wavread(waveFile);
y=y*2^nbits/2;
frameSize=256;
overlap=0;
framedY=buffer(y, frameSize, overlap);
frameNum=size(framedY, 2);
volume=sum(abs(framedY));
voicedIndex=find(volume>frameSize*3);
startIndex=voicedIndex(1);
endIndex=voicedIndex(end);
framedY2=framedY(:, startIndex:endIndex);
PP=setPitchPrm4dp;

allComputed=[];
weight=0:10:100;
for i=weight
	fprintf('%d:\n', i);
	PP.ppcIndexDiffWeight=i;
	[pitch, A, path]=frame2pitchByDp(framedY2, fs, 0, PP);
	computed=zeros(frameNum,1);
	computed(startIndex:endIndex)=pitch(:);
	allComputed=[allComputed, computed];
end
plot(1:frameNum, allComputed);

% 標準答案
load mountain.pitch
desired=mountain;
distance=[];
for i=1:size(allComputed, 2)
	distance(i)=sum(abs(allComputed(:,i)-desired));
end
figure
plot(weight, distance, 'o-');