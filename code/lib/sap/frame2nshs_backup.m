function [nshs, freq]=frame2nshs(frameMat, param, plotOpt)
% frame2nshs: Conversion from a frame to normalized sub-harmonic sum, for pitch tracking
%
%	Usage:
%		[nshs, freq]=frame2nshs(frameMat, param, plotOpt)
%
%	For example:
%		waveFile='soo.wav';
%		startIndex=15000;
%		frameSize=512;
%		[y, param.fs, nbits] = wavread(waveFile);
%		endIndex=startIndex+frameSize-1;
%		frame=y(startIndex:endIndex);
%		param.zeroPaddedFrameSize=16*frameSize;
%		param.polyOrder=10;
%		param.weight=0.8;
%		[nshs, freq]=frame2nshs(frame, param, 1);
%		[maxValue, maxIndex]=max(nshs);
%		fprintf('Pitch freq = %f\n', freq(maxIndex));
%
%	Another example:
%		waveFile='soo.wav';
%		[y, param.fs, nbits] = wavread(waveFile);
%		frameSize=512;
%		overlap=0;
%		frameMat=buffer2(y, frameSize, overlap);
%		[nshs, freq]=frame2nshs(frameMat, [], 1);

%	Roger Jang, 20101017

if nargin<1; selfdemo; return; end
[frameSize, frameNum]=size(frameMat);
if frameSize==1, frameMat=frameMat(:); frameSize=frameNum; end	% Single frame
if nargin<2 | isempty(param)
	param.fs=8000;
	param.zeroPaddedFrameSize=16*length(frameSize);
	param.polyOrder=10;			% Order for polynomial fitting for removing the trend of power spectrum
	param.weight=0.8;			% Weighting factor for summation of power spectrum
	param.pitchRange=[30, 84];
end
if nargin<3; plotOpt=0; end

frameMat2=frameMat;
for i=1:frameNum
	frameMat2(:,i)=frameMat(:,i).*hamming(frameSize);	% Hamming windowing
end

if param.zeroPaddedFrameSize>frameSize
	frameMat2=[frameMat2; zeros(param.zeroPaddedFrameSize-frameSize, frameNum)];		% zero padding
end

[mag, phase, freq, pSpec]=fftOneSide(frameMat2, param.fs, 0);
specLen=length(freq);
pSpec2=pSpec;
if param.polyOrder>=0
	for i=1:frameNum
		pSpecTrend=polyval(polyfit((1:specLen)', pSpec(:,i), param.polyOrder), (1:specLen)');
		pSpec2(:,i)=pSpec(:,i)-pSpecTrend;		% Removing the trend
	end
end

freqRange=pitch2freq(param.pitchRange);
%fprintf('freqRange=%s\n', mat2str(freqRange));
temp=find(freq>=freqRange(1)); index1=temp(1);
temp=find(freq<=freqRange(2)); index2=temp(end);

weightVec=(param.weight).^(0:length(freq)-1)';
nshs=nan*ones(length(freq), frameNum);
for j=1:frameNum
	for i=index1:index2
		step=i-1;
		toBeAdded=pSpec2(i:step:end, j);
		weight=weightVec(1:length(toBeAdded));
		nshs(i, j)=sum(toBeAdded.*weight)/sum(weight);
	end
end

if plotOpt 
	if frameNum==1
		frame=frameMat(:);
		plotNum=4;
		subplot(plotNum,1,1);
		plot(1:length(frame), frame, '.-'); axis tight;
		title('Frame'); xlabel('Samples');

		subplot(plotNum,1,2);
		plot(freq, pSpec, freq, pSpecTrend, 'r');
		set(gca, 'xlim', [-inf inf]);
		title('Power spectrum'); xlabel('Freq (Hz)');

		subplot(plotNum,1,3);
		plot(freq, pSpec2);
		set(gca, 'xlim', [-inf inf]);
		title('Power spectrum'); xlabel('Freq (Hz)');

		subplot(plotNum,1,4);
		plot(freq, nshs); axis tight
		line(freq(index1)*[1 1], [min(nshs), max(nshs)], 'color', 'r');
		line(freq(index2)*[1 1], [min(nshs), max(nshs)], 'color', 'r');
		set(gca, 'xlim', [-inf inf]);
		title('Harmonic product spectrum'); xlabel('Freq (Hz)');
		[maxValue, maxIndex]=max(nshs);
		line(freq(maxIndex), maxValue, 'color', 'r', 'marker', 'o');
	else
		subplot(1,2,1)
		imagesc(mag); axis xy; axis image; colorbar
		title('Magnitude spectrum');
		subplot(1,2,2)
		imagesc(nshs); axis xy; axis image; colorbar
		title('NSHS map');
	end
end

% ====== Self demo
function selfdemo
waveFile='soo.wav';
startIndex=15000;
frameSize=512;
[y, param.fs, nbits] = wavread(waveFile);
endIndex=startIndex+frameSize-1;
frame=y(startIndex:endIndex);
param.zeroPaddedFrameSize=16*frameSize;
param.polyOrder=10;
param.weight=0.8;
param.pitchRange=[30, 84];
[nshs, freq]=feval(mfilename, frame, param, 1);
[maxValue, maxIndex]=max(nshs);
fprintf('Pitch freq = %f\n', freq(maxIndex));
