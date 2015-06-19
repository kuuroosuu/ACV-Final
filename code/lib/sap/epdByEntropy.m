function [epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, entropy, ve] = epdByVolEntropy(wObj, epdPrm, plotOpt)
% epdByVolHod: EPD based on volume and HOD (high-order difference)
%
%	Usage:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, ve] = epdByVolHod(wObj, epdPrm, plotOpt)
%
%	Description:
%		[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec, volume, hod, ve] = epdByVolHod(wObj, epdPrm, plotOpt)
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
%		wObj = waveFile2obj(waveFile);
%		epdPrm=epdPrmSet(wObj.fs);
%		plotOpt = 1;
%		[epInSampleIndex, epInFrameIndex, soundSegment] = epdByEntropy(wObj, epdPrm, plotOpt);

%	Roger Jang, 20110816

if nargin<1, selfdemo; return; end
if nargin<2, fs=16000; end
if nargin<3, nbits=16; end
if nargin<4 | isempty(epdPrm), epdPrm=epdPrmSet(fs); end
if nargin<5, plotOpt=0; end

if size(y, 2)~=1, error('y is not mono!'); end

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
entropy=frame2entropy(frameMat, fs);

% ====== Compute ve thresholds
ve=volume*epdPrm.volWeight+(1-epdPrm.volWeight)*entropy;
% Another method
%volume=volume/max(volume);
entropy=(entropy-min(entropy))/(max(entropy)-min(entropy));
%ve=(volume+(1-entropy))/2;
ve=volume.*(1-entropy);
%ve=volume.*entropy;
temp=sort(ve);
index=round(frameNum*epdPrm.veMinMaxPercentile/100); if index==0, index=1; end
veMin=temp(index);
veMax=temp(frameNum-index+1);			% To avoid unvoiced sounds
veTh=(veMax-veMin)*epdPrm.veRatio+veMin;
%fprintf('veMin=%g, veMax=%g, veTh=%g\n', veMin, veMax, veTh);


epdPrm.fs=fs;
[epInSampleIndex, epInFrameIndex, soundSegment, zeroOneVec]=epdBySingleCurve(ve, epdPrm);

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
	plot(frameTime, volume, '.-', frameTime, entropy, '.-');
	legend('Volume', 'Entropy');
	axis tight;
	ylabel('Volume & Entropy');
	title('Volume & Entropy');

	subplot(3,1,3);
	plot(frameTime, ve, '.-');
	axis tight;
	line([min(frameTime), max(frameTime)], veTh*[1 1], 'color', 'r');
	line([min(frameTime), max(frameTime)], veMin*[1 1], 'color', 'c');
	line([min(frameTime), max(frameTime)], veMax*[1 1], 'color', 'k');
	for i=1:length(soundSegment)
		line(frameTime(soundSegment(i).beginFrame)*[1 1], [0, max(ve)], 'color', 'g');
		line(frameTime(soundSegment(i).endFrame)*[1 1], [0, max(ve)], 'color', 'm');
	end
	ylabel('VE');
	title('VE');
	
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
