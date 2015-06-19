function [epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, volTh] = epdByVol(wObj, epdPrm, plotOpt)
% epdByVol: EPD based on volume only
%
%	Usage:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume] = epdByVol(wObj, epdPrm, plotOpt)
%
%	Description
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume] = epdByVol(wObj, epdPrm, plotOpt)
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
%		out=epdByVol(wObj, epdPrm, plotOpt);

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

y=double(y);					% convert to double data type
frameMat=buffer2(y, frameSize, overlap);	% frame blocking
frameMat=frameZeroMean(frameMat, 2);		% zero justification
frameNum=size(frameMat, 2);					% no. of frames
volume=frame2volume(frameMat);				% compute volume
temp=sort(volume);
index=round(frameNum*epdPrm.vMinMaxPercentile/100); if index==0, index=1; end
volMin=temp(index);
volMax=temp(frameNum-index+1);			% To avoid unvoiced sounds
volTh=(volMax-volMin)*epdPrm.volRatio+volMin;	% compute volume threshold

% ====== Identify voiced part that's larger than volTh2
soundSegment=segmentFind(volume>volTh);

% ====== Delete short sound clips
index = [];
for i=1:length(soundSegment),
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
	subplot(2,1,1);
	time=(1:length(y))/fs;
	frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	plot(time, y);
	for i=1:length(soundSegment)
		line(frameTime(soundSegment(i).beginFrame)*[1 1], 2^nbits/2*[-1, 1], 'color', 'g');
		line(frameTime(soundSegment(i).endFrame)*[1 1], 2^nbits/2*[-1, 1], 'color', 'm');
	end
	axisLimit=[min(time) max(time) -2^nbits/2, 2^nbits/2];
	if -1<=min(y) & max(y)<=1
		axisLimit=[min(time) max(time) -1, 1];
	end
	axis(axisLimit);
	ylabel('Amplitude');
	title('Waveform');

	subplot(2,1,2);
	plot(frameTime, volume, '.-');
	if all(volume)>=0
		axis([-inf inf 0 inf]);
	else
		axis tight;
	end
	line([min(frameTime), max(frameTime)], volTh*[1 1], 'color', 'r');
	line([min(frameTime), max(frameTime)], volMin*[1 1], 'color', 'c');
	line([min(frameTime), max(frameTime)], volMax*[1 1], 'color', 'k');
	for i=1:length(soundSegment)
		line(frameTime(soundSegment(i).beginFrame)*[1 1], [0, max(volume)], 'color', 'g');
		line(frameTime(soundSegment(i).endFrame)*[1 1], [0, max(volume)], 'color', 'm');
	end
	ylabel('Volume');
	title('Volume');
	
	U.y=y; U.fs=fs;
	if max(U.y)>1, U.y=U.y/(2^nbits/2); end
	if ~isempty(epInSampleIndex)
		U.voicedY=U.y(epInSampleIndex(1):epInSampleIndex(end));
	else
		U.voicedY=[];
	end
	set(gcf, 'userData', U);
	uicontrol('string', 'Play all', 'callback', 'U=get(gcf, ''userData''); sound(U.y, U.fs);');
	uicontrol('string', 'Play detected', 'callback', 'U=get(gcf, ''userData''); sound(U.voicedY, U.fs);', 'position', [100, 20, 100, 20]);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
