function onset=wave2onset(wObj, plotOpt)

if nargin<1, selfdemo; return; end

y=wObj.signal;
fs=wObj.fs;
frameSize=256;
overlap=200;
maxTapPerSec=10;
oneSideFrameNum=round(fs/maxTapPerSec/(frameSize-overlap)/2);

opt=wave2volume('defaultOpt');
opt.frameSize=frameSize;
opt.overlap=overlap;
volume=wave2volume(wObj, opt, 1);
frameNum=length(volume);
localIndex=localMax(volume);

localIndexVol=volume(logical(localIndex));
sortedLocalIndexVol=sort(localIndexVol, 'descend');
pos=round(maxTapPerSec*length(y)/fs);
volTh=sortedLocalIndexVol(min(pos, end));
localIndex2=localIndex.*(volume>volTh);

temp=0*(1:length(volume))';
for i=1:frameNum
	if localIndex2(i)<1
		continue;
	end
	center=i;
	startIndex=max(center-oneSideFrameNum, 1);
	endIndex=min(center+oneSideFrameNum, length(volume));
	theVolume=volume(startIndex:endIndex);
	theVolume(logical(1-localIndex(startIndex:endIndex)))=-inf;
	[junk, index]=max(theVolume);
	if i==startIndex+index-1
		temp(i)=1;
	end
end
localIndex3=localIndex(:).*temp;

onset=frame2sampleIndex(find(localIndex3), frameSize, overlap);

if plotOpt
	sampleTime=(1:length(y))/fs;
	frameTime=frame2sampleIndex(1:length(volume), frameSize, overlap)/fs;
	subplot(3,1,1);
	plot(sampleTime, y); set(gca, 'xlim', [-inf inf]);
	axisLimit=axis;
	for i=1:length(onset)
		line(onset(i)*[1, 1]/fs, axisLimit(3:4), 'color', 'r'); 
	end
	xlabel('Time (sec)');
	subplot(3,1,2);
	plot(frameTime, volume, '.-'); set(gca, 'xlim', [-inf inf]);
	localIndex2=find(localIndex2);
	for i=1:length(localIndex2)
		line(frameTime(localIndex2(i)), volume(localIndex2(i)), 'marker', 'o', 'color', 'r');
	end
	line([min(frameTime), max(frameTime)], volTh*[1, 1], 'color', 'k');
	xlabel('Time (sec)');
	subplot(3,1,3);
	plot(frameTime, volume, '.-'); set(gca, 'xlim', [-inf inf]);
	localIndex3=find(localIndex3);
	for i=1:length(localIndex3)
		line(frameTime(localIndex3(i)), volume(localIndex3(i)), 'marker', 'o', 'color', 'r');
	end
	line([min(frameTime), max(frameTime)], volTh*[1, 1], 'color', 'k');
	xlabel('Time (sec)');
end

% ====== Self demo
function selfdemo
waveFile='tapping.wav';
[y, fs, nbits, opts, cueLabel]=wavReadInt(waveFile);
wObj=waveFile2obj(waveFile);
onset=feval(mfilename, wObj, 1);
subplot(3,1,1);
axisLimit=axis;
for i=1:length(cueLabel)
	line(cueLabel(i)*[1]/fs, axisLimit(3), 'color', 'k', 'marker', '^'); 
end