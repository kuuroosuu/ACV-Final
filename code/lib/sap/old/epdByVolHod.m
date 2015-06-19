function [epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, vh] = epdByVolHod(y, fs, nbits, epdPrm, plotOpt)
% epdByVolHod: EPD based on volume and HOD (high-order difference)
%
%	Usage:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, vh] = epdByVolHod(y, fs, nbits, epdPrm, plotOpt)
%
%	Description:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, vh] = epdByVolHod(y, fs, nbits, epdPrm, plotOpt)
%			epInSampleIndex: two-element end-points in sample index
%			epInFrameIndex: two-element end-points in frame index
%			soundSegment: segment of voice activity
%			y: input audio signals
%			fs: sampling rate
%			nbits: no. of bits
%			epdPrm: parameters for EPD
%			plotOpt: 0 for silence operation, 1 for plotting
%
%	Example:
%		waveFile='SingaporeIsAFinePlace.wav';
%		[y, fs, nbits] = wavReadInt(waveFile);
%		epdPrm=epdPrmSet(fs);
%		plotOpt = 1;
%		[epInSampleIndex, epInFrameIndex, soundSegment] = epdByVolHod(y, fs, nbits, epdPrm, plotOpt);

%	Roger Jang, 20070323

if nargin<1, selfdemo; return; end
if nargin<2, fs=16000; end
if nargin<3, nbits=16; end
if nargin<4 | isempty(epdPrm), epdPrm=epdPrmSet(fs); end
if nargin<5, plotOpt=0; end

if size(y, 2)~=1, error('y is not mono!'); end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
minSegment=round(epdPrm.minSegment*fs/(frameSize-overlap));
maxSilBetweenWord=round(epdPrm.maxSilBetweenWord*fs/(frameSize-overlap));
%minLastWordDuration=round(epdPrm.minLastWordDuration*fs/(frameSize-overlap));

% ====== Compute volume/hod
frameMat=buffer2(y, frameSize, overlap);	% frame blocking
frameMat=frameZeroMean(frameMat, 2);
frameNum=size(frameMat, 2);			% no. of frames
volume=frame2volume(frameMat);
hod=frame2ashod(frameMat, epdPrm.diffOrder);

% ====== Compute vh and its thresholds
volume=volume/max(volume); hod=hod/max(hod);	% Normalization before mixing
vh=volume*epdPrm.volWeight+(1-epdPrm.volWeight)*hod;
temp=sort(vh);
index=round(frameNum*epdPrm.vhMinMaxPercentile/100); if index==0, index=1; end
vhMin=temp(index);
vhMax=temp(frameNum-index+1);			% To avoid unvoiced sounds
vhTh=(vhMax-vhMin)*epdPrm.vhRatio+vhMin;
%fprintf('vhMin=%g, vhMax=%g, vhTh=%g\n', vhMin, vhMax, vhTh);

% ====== Identify voiced part that's larger than volTh2
soundSegment=segmentFind(vh>vhTh);

% ====== Delete short sound clips
index = [];
for i=1:length(soundSegment),
	if soundSegment(i).duration<=minSegment
		index = [index, i];
	end
end
soundSegment(index) = [];

% ====== If the sil between the last two segment is too big, delete the last segment
if length(soundSegment)>1
	if soundSegment(end).begin-soundSegment(end-1).end>maxSilBetweenWord
%		if soundSegment(end).duration<=minLastWordDuration		% Not used!
			soundSegment(end)=[];
%		end
	end
end

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

zeroOneVec=0*vh;
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
		soundSegment(i).endSample   = min(length(y), frame2sampleIndex(soundSegment(i).end, frameSize, overlap));
		soundSegment(i).beginFrame = soundSegment(i).begin;
		soundSegment(i).endFrame = soundSegment(i).end;
	end
	soundSegment=rmfield(soundSegment, 'begin');
	soundSegment=rmfield(soundSegment, 'end');
	soundSegment=rmfield(soundSegment, 'duration');
end

% Plotting
if plotOpt,
	subplot(3,1,1);
	time=(1:length(y))/fs;
	frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	plot(time, y);
	for i=1:length(soundSegment)
		line(frameTime(soundSegment(i).beginFrame)*[1 1], 2^nbits/2*[-1, 1], 'color', 'm');
		line(frameTime(soundSegment(i).endFrame)*[1 1], 2^nbits/2*[-1, 1], 'color', 'g');
	end
	axisLimit=[min(time) max(time) -2^nbits/2, 2^nbits/2];
	if -1<=min(y) & max(y)<=1
		axisLimit=[min(time) max(time) -1, 1];
	end
	axis(axisLimit);
	ylabel('Amplitude');
	title('Waveform');
	
	subplot(3,1,2);
	plot(frameTime, volume, '.-', frameTime, hod, '.-');
	legend('Volume', 'HOD');
	axis tight;
	ylabel('Volume & HOD');
	title('Volume & HOD');

	subplot(3,1,3);
	plot(frameTime, vh, '.-');
	axis tight;
	line([min(frameTime), max(frameTime)], vhTh*[1 1], 'color', 'r');
	line([min(frameTime), max(frameTime)], vhMin*[1 1], 'color', 'c');
	line([min(frameTime), max(frameTime)], vhMax*[1 1], 'color', 'k');
	for i=1:length(soundSegment)
		line(frameTime(soundSegment(i).beginFrame)*[1 1], [0, max(vh)], 'color', 'm');
		line(frameTime(soundSegment(i).endFrame)*[1 1], [0, max(vh)], 'color', 'g');
	end
	ylabel('VH');
	title('VH');
	
	U.y=double(y); U.fs=fs;
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
