function [endPoint, volume, zcr, soundSegment] = epdByVolZcr(y, fs, nbits, epdPrm, plotOpt)
% epdByVolZcr: End point detection based on volume and zero-crossing rate
%	Usage: [endPoint, volume, zcr, soundSegment] = epdByVolZcr(y, fs, nbits, epdPrm, plotOpt)

%	Roger Jang, 20041118

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, nbits=8; end
if nargin<4 | isempty(epdPrm)
	epdPrm.frameSize = round(fs/31.25);		% fs=8000 ===> frameSize=256;
	epdPrm.overlap = round(epdPrm.frameSize/2);
	epdPrm.volRatio1=0.1;		% lower threshold
	epdPrm.volRatio2=0.2;		% upper threshold
	epdPrm.zcrRatio=0.1;
	epdPrm.zcrShiftGain=4;
end
if nargin<5, plotOpt=0; end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
zcrRatio=epdPrm.zcrRatio;

if max(y)-min(y)<=2
	error('The range of the given y vector is too small! Perhaps it should be converted to integer format...');
end

% ====== Zero adjusted
y = y-round(mean(y));

% ====== Take frames
frameMat  = buffer2(y, frameSize, overlap);
frameNum = size(frameMat, 2);			% Number of frames
frameTime = frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;

% ====== Compute volume/zcr
volume=frame2volume(frameMat);
[minVol, index]=min(volume);
shiftAmount=epdPrm.zcrShiftGain*max(abs(frameMat(:,index)));		% shiftAmount is equal to epdPrm.zcrShiftGain times the max. abs. sample within the frame of min. volume
shiftAmount=max(shiftAmount, 2);	
zcr=frame2zcr(frameMat, 1, shiftAmount);

% ====== Compute volume/zcr thresholds
maxVol=max(volume);
volTh2=(maxVol-minVol)*epdPrm.volRatio2+minVol;
volTh1=(maxVol-minVol)*epdPrm.volRatio1+minVol;
zcrTh=max(zcr)*zcrRatio;

% ====== Identify voiced part that's larger than volTh2
voicedIndex = volume>=volTh2;
soundSegment=findSegment(voicedIndex);

% ====== Delete short sound clips
%index = [];
%for i=1:length(soundSegment),
%	if (soundSegment(i).end-soundSegment(i).begin)<4
%		index = [index, i];
%	end
%end
%soundSegment(index) = [];

% ====== Expansion 1: Expand end points to volume level1 (lower level)
for i=1:length(soundSegment),
	head = soundSegment(i).begin;
	while (head-1)>=1 & volume(head-1)>=volTh1,
		head=head-1;
	end
	soundSegment(i).begin = head;
	tail = soundSegment(i).end;
	while (tail+1)<=length(volume) & volume(tail+1)>=volTh1,
		tail=tail+1;
	end
	soundSegment(i).end = tail;
end
% ====== Expansion 2: Expand end points to include high zcr region
for i=1:length(soundSegment),
	head = soundSegment(i).begin;
	while (head-1)>=1 & zcr(head-1)>zcrTh			% Extend at beginning
		head=head-1;
	end
	soundSegment(i).begin = head;
	tail = soundSegment(i).end;
	while (tail+1)<=length(zcr) & zcr(tail+1)>zcrTh		% Extend at ending
		tail=tail+1;
	end
	soundSegment(i).end = tail;
end

% ====== Delete repeated sound segments
index = [];
for i=1:length(soundSegment)-1,
	if soundSegment(i).begin==soundSegment(i+1).begin & soundSegment(i).end==soundSegment(i+1).end,
		index=[index, i];
	end
end
soundSegment(index) = [];

% ====== Transform sample-point-based index
if length(soundSegment) ~=0,
	for i=1:length(soundSegment),
		soundSegment(i).beginSample = frame2sampleIndex(soundSegment(i).begin, frameSize, overlap);
		soundSegment(i).endSample   = min(length(y), frame2sampleIndex(soundSegment(i).end, frameSize, overlap));
		soundSegment(i).beginFrame = soundSegment(i).begin;
		soundSegment(i).endFrame = soundSegment(i).end;
	end
	endPoint=[soundSegment(1).beginSample, soundSegment(end).endSample];	% 取頭尾
else
	endPoint = [];
end
soundSegment=rmfield(soundSegment, 'begin');
soundSegment=rmfield(soundSegment, 'end');
soundSegment=rmfield(soundSegment, 'duration');

% ====== 如果相鄰的 segment 間隔不到0.2秒，則合併
gap=[];
for i=1:length(soundSegment)-1
	gap(i)=(soundSegment(i+1).beginSample-soundSegment(i).endSample)/fs;
end
group=findSegment(gap<0.2);
deletedIndex=[];
for i=1:length(group)
	soundSegment(group(i).begin).endSample=soundSegment(group(i).end+1).endSample;
	soundSegment(group(i).begin).endFrame=soundSegment(group(i).end+1).endFrame;
	deletedIndex=[deletedIndex, (group(i).begin+1):(group(i).end+1)];
end
soundSegment(deletedIndex)=[];

if plotOpt,
	axes1H=subplot(4,1,1);
	time=(1:length(y))/fs;
	plot(time, y);
	axis([min(time), max(time), -2^nbits/2, 2^nbits/2]);
	ylabel('Amplitude'); title('Waveform'); grid on
	% Plot end points
	yBound=[-2^nbits/2, 2^nbits/2];
	for i=1:length(soundSegment),
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes2H=subplot(4,1,2);
	frameTime = frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	plot(frameTime, volume, '.-');
	line([min(frameTime), max(frameTime)], volTh1*[1 1], 'color', 'r');
	line([min(frameTime), max(frameTime)], volTh2*[1 1], 'color', 'r');
	axis tight
	ylabel('Volume'); title('Volume'); grid on
	% Plot end points
	yBound = [min(volume) max(volume)];
	for i=1:length(soundSegment),
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes3H=subplot(4,1,3);
	plot(frameTime, zcr, '.-');
	line([min(frameTime), max(frameTime)], zcrTh*[1 1], 'color', 'c');
	axis([min(frameTime), max(frameTime), 0, max(zcr)]);
	ylabel('ZCR'); title('Zero crossing rate'); grid on
	% Plot end points
	yBound = [0 max(zcr)];
	for i=1:length(soundSegment),
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes4H=subplot(4,1,4);
	voicedIndex=endPoint(1):endPoint(2);
	voicedTime=time(voicedIndex);
	voicedY=y(voicedIndex);
	voicedH=plot(voicedTime, voicedY);
	axis([time(endPoint(1)), time(endPoint(2)), -2^nbits/2, 2^nbits/2]);
	ylabel('Amplitude'); title('Voiced waveform'); grid on
	
	U.y=y; U.fs=fs; U.nbits=nbits;
	U.axes1H=axes1H; U.axes2H=axes2H; U.axes3H=axes3H; U.axes4H=axes4H;
	U.voicedIndex=voicedIndex; U.voicedH=voicedH;
	U.voicedY=voicedY; U.voicedTime=voicedTime;
	set(gcf, 'userData', U);
	uicontrol('string', 'Play all', 'callback', 'U=get(gcf, ''userData''); sound(U.y/(2^U.nbits/2), U.fs);');
	uicontrol('string', 'Play voiced', 'callback', 'U=get(gcf, ''userData''); sound(U.voicedY/(2^U.nbits/2), U.fs);', 'position', [100, 20, 60, 20]);

	% Play the segmented sound
%	head = soundSegment(1).beginFrame*(frameSize-overlap);
%	tail = min(length(y), soundSegment(end).end*(frameSize-overlap));
%	thisY = y(head:tail);
%	fprintf('His return to hear the cutted sound %g:', i);
%	pause;
%	fprintf('\n');
%	wavplay(thisY, fs, 'sync');
%	fprintf('\n');
end

% ====== Self demo
function selfdemo
waveFile='SingaporeIsAFinePlace.wav';
plotOpt = 1;
[y, fs, nbits] = wavReadInt(waveFile);
endPoint = feval(mfilename, y, fs, nbits, [], plotOpt);