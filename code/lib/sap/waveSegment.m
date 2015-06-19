function cutPointSampleIndex=waveSegment(y, fs, minDuration, plotOpt)
% waveSegment: Segmentation of wave using its volume
%	Usage: cutPointSampleIndex=waveSegment(y, fs, minDuration, plotOpt)

%	Roger Jang

if nargin<1, selfdemo; return; end
if nargin<2, fs=16000; end
if nargin<3, minDuration=length(y)/fs/10; end
if nargin<4, plotOpt=0; end

frameSize=fs/8000*512;
overlap=frameSize/2;
frameMat=buffer2(y, frameSize, overlap);
volume=frame2volume(frameMat);
minLen=round(minDuration*fs/(frameSize-overlap));
cutPoint=volumeSegment(volume, minLen);
cutPointSampleIndex=frame2sampleIndex(cutPoint, frameSize, overlap);

if plotOpt
	subplot(2,1,1);
	plot((1:length(y))/fs, y);
	set(gca, 'xlim', inf*[-1, 1]);
	axisLimit=axis;
	for i=1:length(cutPointSampleIndex)
		line(cutPointSampleIndex(i)*[1 1]/fs, axisLimit(3:4), 'color', 'r');
	end
	
	subplot(2,1,2);
	plot(volume); axis tight
	axisLimit=axis;
	for i=1:length(cutPoint)
		line(cutPoint(i)*[1 1], axisLimit(3:4), 'color', 'r');
	end
	xlabel('Frame index');
	ylabel('Volume');
end

% ====== Selfdemo
function selfdemo
waveFile='query0001.wav';
[y, fs, nbits]=wavread(waveFile);
cutPointSampleIndex=feval(mfilename, y, fs, 2, 1);
% Playback of segments
temp=[1, cutPointSampleIndex, length(y)];
segmentNum=length(temp)-1;
segmentDuration=diff(temp)/fs;
fprintf('Duration of segments = %s\n', mat2str(segmentDuration));
for i=1:segmentNum
	thisWave=y(temp(i):temp(i+1));
	fprintf('Press any key to play segment %d/%d (%f sec)...', i, segmentNum, segmentDuration(i)); pause; fprintf('\n');
	sound(thisWave, fs);
end