function [endPoint, volume, zcr, soundSegment] = endPointDetect03(y, fs, nbits, plotOpt, epdPrm)
% endPointDetect: End point detection based on volume only, for PMP recorded sound
%	Usage: [endPoint, volume, zcr, soundSegment] = endPointDetect03(y, fs, nbits, plotOpt, epdPrm)

%	Roger Jang, 20041118

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, nbits=8; end
if nargin<4, plotOpt=0; end
if nargin<5,
	epdPrm.frameSize = round(fs*480/16000);	% fs=16000 ===> frameSize=480;
	epdPrm.overlap = round(epdPrm.frameSize/3*2);
	epdPrm.zcrRatio = 0.25;
	epdPrm.volRatio1=10;
	epdPrm.volRatio2=20;
end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
zcrRatio=epdPrm.zcrRatio;

if max(y)-min(y)<=2
	error('The range of the given y vector is too small! Perhaps it should be converted to integer format...');
end

% ====== Zero adjusted
%y = y-round(mean(y));

% ====== Take frames
framedY  = buffer2(y, frameSize, overlap);
frameNum = size(framedY, 2);			% Number of frames
frameTime = frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;

% ====== Compute volume/zcr
volume=frame2volume(framedY);
[minVol, index]=min(volume);
shiftAmount=2*max(abs(framedY(:,index)));	% 以最小音量之音框中的訊號最大絕對值的兩倍為平移量來進行ZCR
shiiftAmount=max(shiftAmount, 2);	
zcr=frame2zcr(framedY, shiftAmount);

% ====== Compute volume/zcr thresholds
volTh1=epdPrm.volRatio1*minVol;
volTh2=epdPrm.volRatio2*minVol;
maxVol=max(volume);
if volTh2>max(volume)/2			% 避免雜訊太大，導致 volTh2 大於音量最大值的錯誤情況
	fprintf('volTh2 too big!\n');
	volTh1=maxVol/5;
	volTh2=2*volTh1;
end
if maxVol>128*volTh2			% 避免雜訊太小
	volTh1=maxVol/128;
	volTh2=2*volTh1;
end
zcrTh = max(zcr)*zcrRatio;

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

% ====== Expand end points to volume level1 (lower level)
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

if plotOpt,
	axes1H=subplot(4,1,1);
	time=(1:length(y))/fs;
	plot(time, y);
	axis([min(time), max(time), -2^nbits/2, 2^nbits/2]);
	ylabel('Amplitude');
	title('Waveform');
	% Plot end points
	yBound=[-2^nbits/2, 2^nbits/2];
	for i=1:length(soundSegment),
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes2H=subplot(4,1,2);
	frameTime = frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	plot(frameTime, volume, '.-');
	line([min(frameTime), max(frameTime)], volTh1*[1 1], 'color', 'c');
	line([min(frameTime), max(frameTime)], volTh2*[1 1], 'color', 'c');
	axis tight
	ylabel('Volume');
	title('Volume');
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
	ylabel('ZCR');
	title('Zero crossing rate');
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
	ylabel('Amplitude');
	title('Voiced waveform');
	
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
waveFile='清華大學資訊系.wav';
waveFile='yankee doodle_不詳_0.wav';
plotOpt = 1;
[y, fs, nbits] = wavReadInt(waveFile);
endPoint = feval(mfilename, y, fs, nbits, plotOpt);