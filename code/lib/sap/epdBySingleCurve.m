function [epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec]=epdBySingleCurve(epdCurve, epdPrm, plotOpt)
% epdBySingleCurve: EPD based on a single curve (volume, HOD, etc)
%
%	Usage:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec] = epdBySingleCurve(epdCurve, epdPrm, plotOpt)
%
%	Description:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec] = epdBySingleCurve(epdCurve, epdPrm, plotOpt)
%			epInSampleIndex: two-element end-points in sample index
%			epInFrameIndex: two-element end-points in frame index
%			soundSegment: segment of voice activity
%			epdCurve: curve for EPD
%			epdPrm: parameters for EPD
%			plotOpt: 0 for silence operation, 1 for plotting
%
%	Example:
%		waveFile='SingaporeIsAFinePlace.wav';
%		[y, fs, nbits] = wavread(waveFile);
%		epdPrm=epdPrmSet(fs);
%		epdPrm.fs=fs;
%		plotOpt=1;
%		frameMat=enframe(y, epdPrm.frameSize, epdPrm.overlap);	% frame blocking
%		frameMat=frameZeroMean(frameMat, 2);
%		vol=frame2volume(frameMat);
%		hod=frame2ashod(frameMat, epdPrm.diffOrder);
%		vol=vol/max(vol); hod=hod/max(hod);	% Normalization before mixing
%		vh=vol*epdPrm.volWeight+(1-epdPrm.volWeight)*hod;
%		[epInSampleIndex, epInFrameIndex, soundSegment]=epdBySingleCurve(vh, epdPrm, plotOpt);

%	Roger Jang, 20110828

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
fs=epdPrm.fs;
% Convert into the unit of frames
minSegment=round(epdPrm.minSegment*fs/(frameSize-overlap));
maxSilBetweenSegment=round(epdPrm.maxSilBetweenSegment*fs/(frameSize-overlap));
%minLastWordDuration=round(epdPrm.minLastWordDuration*fs/(frameSize-overlap));
frameNum=length(epdCurve);

% ====== Compute vh and its thresholds
temp=sort(epdCurve);
index=round(frameNum*epdPrm.vhMinMaxPercentile/100); if index==0, index=1; end
%fprintf('index=%d, epdPrm.vhThPercentile=%f\n', index, epdPrm.vhThPercentile);
epdCurveMin=temp(index);
epdCurveMax=temp(frameNum-index+1);			% To avoid unvoiced sounds
epdCurveTh=(epdCurveMax-epdCurveMin)*epdPrm.vhRatio+epdCurveMin;
%fprintf('epdCurveMin=%g, epdCurveMax=%g, epdCurveTh=%g\n', epdCurveMin, epdCurveMax, epdCurveTh);

% ====== Identify voiced part that's larger than volTh2
soundSegment=segmentFind(epdCurve>epdCurveTh);

% ====== Delete short sound clips
index = [];
for i=1:length(soundSegment),
	if soundSegment(i).duration<=minSegment
		index = [index, i];
	end
end
soundSegment(index) = [];

% ====== Assign values to segment
for i=1:length(soundSegment)
	soundSegment(i).value=epdCurve(soundSegment(i).begin:soundSegment(i).end);
end

% Find the max gap between segment...
if length(soundSegment)>=2
	bb=[soundSegment.begin];
	ee=[soundSegment.end];
	gap=bb(2:end)-ee(1:end-1);
	[maxGap, maxIndex]=max(gap);
	if maxGap>maxSilBetweenSegment
	%	weight1=sum([soundSegment(1:maxIndex).duration]);
	%	weight2=sum([soundSegment(maxIndex+1:end).duration]);
		weight1=sum([soundSegment(1:maxIndex).value]);
		weight2=sum([soundSegment(maxIndex+1:end).value]);
		if weight1>weight2
			soundSegment(maxIndex+1:end)=[];
		else
			soundSegment(1:maxIndex)=[];
		end
	end
end

%{
% ====== If the sil between the last two segment is too big and if the last segment is too short, delete the last segment
if length(soundSegment)>=2
	if soundSegment(end).begin-soundSegment(end-1).end>maxSilBetweenSegment
	%	if soundSegment(end).duration<=2*minSegment
			soundSegment(end)=[];
	%	end
	end
end
% ====== If the sil between the first two segment is too big and if the first segment is too short, delete the last segment
if length(soundSegment)>=2
	if soundSegment(2).begin-soundSegment(1).end>maxSilBetweenSegment
	%	if soundSegment(1).duration<=2*minSegment
			soundSegment(1)=[];
	%	end
	end
end
%}

% ====== Check if segment 1 is noise
%if length(soundSegment)>=2
%	sil=soundSegment(2).begin-soundSegment(1).end-1;
%	if (sil>epdPrm.minSil) && (2*soundSegment(1).duration<soundSegment(2).duration)
%		soundSegment(1)=[];
%	end
%end
% ====== Check if the last segment is noise
%if length(soundSegment)>=2
%	sil=soundSegment(end).begin-soundSegment(end-1).end-1;
%	if (sil>epdPrm.minSil) && (soundSegment(end-1).duration>2*soundSegment(end).duration)
%		soundSegment(end)=[];
%	end
%end

% Use segment2 for further processing. Keep segment for plotting
soundSegment2=soundSegment;
% If a sil is longer than its neighboring segment, delete the segment
while 0
	silDuration=[];
	for i=1:length(soundSegment2)-1
		silDuration(i)=soundSegment2(i+1).begin-soundSegment2(i).end-1;
	end
	if isempty(silDuration), break; end
	[maxSil, index]=max(silDuration);
%	if maxSil>epdPrm.minSil	% max silence too long ===> delete one of its neighboring segment
	if maxSil>min(soundSegment2(index).duration, soundSegment2(index+1).duration)		% max silence too long ===> delete one of its neighboring segment
		if soundSegment2(index).duration<soundSegment2(index+1).duration
			soundSegment2(index)=[];
		else
			soundSegment2(index+1)=[];
		end
	else
		break;
	end
end

% Entending endpoints
%soundSegment2(1).begin=max(soundSegment2(1).begin-epdPrm.extendNum, 1);
%soundSegment2(end).end=min(soundSegment2(end).end+epdPrm.extendNum, frameNum);

zeroOneVec=logical(0*epdCurve);
for i=1:length(soundSegment2)
	for j=soundSegment2(i).begin:soundSegment2(i).end
		zeroOneVec(j)=1;
	end
end

if isempty(soundSegment2)
	epInSampleIndex=[];
	epInFrameIndex=[];
	fprintf('Warning: No segment found in %s.m.\n', mfilename);
else
	epInFrameIndex=[soundSegment2(1).begin, soundSegment2(end).end];
	epInSampleIndex=frame2sampleIndex(epInFrameIndex, frameSize, overlap);		% conversion from frame index to sample index
	for i=1:length(soundSegment),
		soundSegment(i).beginSample = frame2sampleIndex(soundSegment(i).begin, frameSize, overlap);
		soundSegment(i).endSample   = frame2sampleIndex(soundSegment(i).end, frameSize, overlap);
		soundSegment(i).beginFrame = soundSegment(i).begin;
		soundSegment(i).endFrame = soundSegment(i).end;
	end
	soundSegment=rmfield(soundSegment, 'begin');
	soundSegment=rmfield(soundSegment, 'end');
%	soundSegment=rmfield(soundSegment, 'duration');
end

% Plotting
if plotOpt
	frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	plot(frameTime, epdCurve, '.-');
	axis tight;
	line([min(frameTime), max(frameTime)], epdCurveTh*[1 1], 'color', 'r');
	line([min(frameTime), max(frameTime)], epdCurveMin*[1 1], 'color', 'c');
	line([min(frameTime), max(frameTime)], epdCurveMax*[1 1], 'color', 'k');
	for i=1:length(soundSegment)
		line(frameTime(soundSegment(i).beginFrame)*[1 1], [0, max(epdCurve)], 'color', 'g');
		line(frameTime(soundSegment(i).endFrame)*[1 1], [0, max(epdCurve)], 'color', 'm');
	end
	ylabel('epdCurve');
	title('epdCurve');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
