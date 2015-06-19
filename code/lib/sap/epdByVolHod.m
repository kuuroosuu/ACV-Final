function [epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, vh] = epdByVolHod(wObj, epdPrm, plotOpt)
% epdByVolHod: EPD based on volume and HOD (high-order difference)
%
%	Usage:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, vh] = epdByVolHod(wObj, epdPrm, plotOpt)
%
%	Description:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, vh] = epdByVolHod(wObj, epdPrm, plotOpt)
%			epInSampleIndex: two-element end-points in sample index
%			epInFrameIndex: two-element end-points in frame index
%			soundSegment: segment of voice activity
%			zeroOneVec: zero-one vector for each frame
%			volume: volume
%			wObj: wave object
%			epdPrm: parameters for EPD
%			plotOpt: 0 for silence operation, 1 for plotting
%
%	Example:
%		waveFile='SingaporeIsAFinePlace.wav';
%		wObj = waveFile2obj(waveFile);
%		epdPrm=epdPrmSet(wObj.fs);
%		plotOpt = 1;
%		[epInSampleIndex, epInFrameIndex, soundSegment] = epdByVolHod(wObj, epdPrm, plotOpt);

%	Roger Jang, 20070323

if nargin<1, selfdemo; return; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name
if nargin<2 | isempty(epdPrm), epdPrm=epdPrmSet(fs); end
if nargin<3, plotOpt=0; end

y=wObj.signal; fs=wObj.fs; nbits=wObj.nbits;
if size(y, 2)~=1, error('Wave is not mono!'); end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
minSegment=round(epdPrm.minSegment*fs/(frameSize-overlap));
maxSilBetweenSegment=round(epdPrm.maxSilBetweenSegment*fs/(frameSize-overlap));
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

epdPrm.fs=fs;
[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec]=epdBySingleCurve(vh, epdPrm);

% ====== Plotting
if plotOpt,
	subplot(3,1,1);
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
		line(frameTime(soundSegment(i).beginFrame)*[1 1], [0, max(vh)], 'color', 'g');
		line(frameTime(soundSegment(i).endFrame)*[1 1], [0, max(vh)], 'color', 'm');
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
