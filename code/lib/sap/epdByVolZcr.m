function [epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume] = epdByVolZcr(wObj, epdPrm, plotOpt)
% epdByVolZcr: EPD based on volume and zero-crossing rate (ZCR)
%
%	Usage:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume] = epdByVolZcr(wObj, epdPrm, plotOpt)
%
%	Description:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume] = epdByVolZcr(wObj, epdPrm, plotOpt)
%			epInSampleIndex: two-element end-points in sample index
%			epInFrameIndex: two-element end-points in frame index
%			soundSegment: resulting sound segments
%			zeroOneVec: zero-one vector for each frame
%			volume: volume
%			wObj: wave object
%			epdPrm: parameters for EPD
%			plotOpt: 0 for silence operation, 1 for plotting
%
%	Example:
%		waveFile='SingaporeIsAFinePlace.wav';
%		wObj=waveFile2obj(waveFile);
%		epdPrm=epdPrmSet(wObj.fs);
%		plotOpt=1;
%		out=epdByVolZcr(wObj, epdPrm, plotOpt);

%	Roger Jang, 20040413, 20070320

if nargin<1, selfdemo; return; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name
if nargin<2 | isempty(epdPrm), epdPrm=epdPrmSet(fs); end
if nargin<3, plotOpt=0; end

y=wObj.signal; fs=wObj.fs; nbits=wObj.nbits;
if size(y, 2)~=1, error('Wave is not mono!'); end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
minSegment=round(epdPrm.minSegment*fs/(frameSize-overlap));		% In terms of no. of frames
maxSilBetweenSegment=round(epdPrm.maxSilBetweenSegment*fs/(frameSize-overlap));		% In terms of no. of frames
%minLastWordDuration=round(epdPrm.minLastWordDuration*fs/(frameSize-overlap));

y = double(y);					% convert to double data type
frameMat=enframe(y, frameSize, overlap);	% frame blocking
frameMat=frameZeroMean(frameMat, 2);
frameNum=size(frameMat, 2);			% no. of frames
volume=frame2volume(frameMat);		% compute volume
temp=sort(volume);
index=round(frameNum*epdPrm.vMinMaxPercentile/100); if index==0, index=1; end
volMin=temp(index);
volMax=temp(frameNum-index+1);			% To avoid unvoiced sounds
volTh1=(volMax-volMin)*epdPrm.volRatio+volMin;	% compute volume threshold
volTh2=(volMax-volMin)*epdPrm.volRatio2+volMin;	% compute volume threshold

% ====== Identify voiced part that's larger than volTh2
soundSegment=segmentFind(volume>volTh1);

% ====== Compute ZCR
[minVol, index]=min(volume);
shiftAmount=epdPrm.zcrShiftGain*max(abs(frameMat(:,index)));		% shiftAmount is equal to epdPrm.zcrShiftGain times the max. abs. sample within the frame of min. volume
%shiftAmount=max(shiftAmount, 2);
shiftAmount=max(shiftAmount, max(frameMat(:))/100);
zcr=frame2zcr(frameMat, 1, shiftAmount);
zcrTh=max(zcr)*epdPrm.zcrRatio;

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

% ====== Delete short sound clips
index = [];
for i=1:length(soundSegment)
	soundSegment(i).duration=soundSegment(i).end-soundSegment(i).begin+1;	% This is necessary since the duration is changed due to expansion
	if soundSegment(i).duration<=minSegment
		index = [index, i];
	end
end
soundSegment(index) = [];

zeroOneVec=0*volume;
for i=1:length(soundSegment)
	for j=soundSegment(i).begin:soundSegment(i).end
		zeroOneVec(j)=1;
	end
end

if isempty(soundSegment)
	epInSampleIndex=[];
	epInFrameIndex=[];
	fprintf('Warning: No sound segment found in %s.m.\n', mfilename);
else
	epInFrameIndex=[soundSegment(1).begin, soundSegment(end).end];
	epInSampleIndex=frame2sampleIndex(epInFrameIndex, frameSize, overlap);		% conversion from frame index to sample index
	for i=1:length(soundSegment),
		soundSegment(i).beginSample = frame2sampleIndex(soundSegment(i).begin, frameSize, overlap);
		soundSegment(i).endSample   = min(length(y), frame2sampleIndex(soundSegment(i).end, frameSize, overlap));
		soundSegment(i).beginFrame = soundSegment(i).begin;
		soundSegment(i).endFrame = soundSegment(i).end;
	end
	soundSegment=rmfield(soundSegment, 'begin');
	soundSegment=rmfield(soundSegment, 'end');
	soundSegment=rmfield(soundSegment, 'duration');
end

% ====== Plotting...
if plotOpt,
	axes1H=subplot(4,1,1);
	time=(1:length(y))/fs;
	plot(time, y);
	axisLimit=[min(time) max(time) -2^nbits/2, 2^nbits/2];
	if -1<=min(y) & max(y)<=1
		axisLimit=[min(time) max(time) -1, 1];
	end
	axis(axisLimit);
	ylabel('Amplitude'); title('Waveform'); grid on
	% Plot end points
	yBound=axisLimit(3:4);
	for i=1:length(soundSegment),
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
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
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
	end

	axes3H=subplot(4,1,3);
	plot(frameTime, zcr, '.-');
	line([min(frameTime), max(frameTime)], zcrTh*[1 1], 'color', 'c');
	axis([min(frameTime), max(frameTime), 0, max(zcr)]);
	ylabel('ZCR'); title('Zero crossing rate'); grid on
	% Plot end points
	yBound = [0 max(zcr)];
	for i=1:length(soundSegment),
		line(frame2sampleIndex(soundSegment(i).beginFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
		line(frame2sampleIndex(  soundSegment(i).endFrame, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
	end

	axes4H=subplot(4,1,4);
	voicedIndex=epInSampleIndex(1):epInSampleIndex(2);
	voicedTime=time(voicedIndex);
	voicedY=y(voicedIndex);
	voicedH=plot(voicedTime, voicedY);
	axis([time(epInSampleIndex(1)), time(epInSampleIndex(2)), -1, 1]);
	ylabel('Amplitude'); title('Voiced waveform'); grid on
	
	U.y=y; U.fs=fs; U.nbits=nbits;
	U.axes1H=axes1H; U.axes2H=axes2H; U.axes3H=axes3H; U.axes4H=axes4H;
	U.voicedIndex=voicedIndex; U.voicedH=voicedH;
	U.voicedY=voicedY; U.voicedTime=voicedTime;
	set(gcf, 'userData', U);
	uicontrol('string', 'Play all', 'callback', 'U=get(gcf, ''userData''); sound(U.y, U.fs);');
	uicontrol('string', 'Play detected', 'callback', 'U=get(gcf, ''userData''); sound(U.voicedY, U.fs);', 'position', [100, 20, 100, 20]);

	% Play the segmented sound
%	head = soundSegment(1).beginFrame*(frameSize-overlap);
%	tail = min(length(y), soundSegment(end).endFrame*(frameSize-overlap));
%	thisY = y(head:tail);
%	fprintf('His return to hear the cutted sound %g:', i);
%	pause;
%	fprintf('\n');
%	sound(thisY, fs);
%	fprintf('\n');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
