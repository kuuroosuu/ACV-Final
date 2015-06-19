function [nshs, specMat2, specTrend]=spec2nshs(specMat, opt, plotOpt)
% spec2nshs: Conversion from a spectrum to normalized sub-harmonic sum
%
%	Usage:
%		[nshs, specMat2, specTrend]=spec2nshs(specMat, opt, plotOpt)
%
%	For example, the following example compute the pitch of NSHS from a single frame:
%		waveFile='soo.wav';
%		startIndex=15000;
%		frameSize=512;
%		[y, opt.fs, nbits] = wavread(waveFile);
%		endIndex=startIndex+frameSize-1;
%		frame=y(startIndex:endIndex);
%		opt.zeroPaddedFrameSize=16*frameSize;
%		opt.polyOrder=10;
%		opt.weight=0.8;
%		[nshs, freq]=frame2nshs(frame, opt, 1);
%		[maxValue, maxIndex]=max(nshs);
%		fprintf('Pitch freq = %f\n', freq(maxIndex));
%
%	Another example:
%		waveFile='soo.wav';
%		[y, opt.fs, nbits] = wavread(waveFile);
%		frameSize=512;
%		overlap=0;
%		frameMat=buffer2(y, frameSize, overlap);
%		[magSpec, phaseSpec, opt.freq, powerSpec]=fftOneSide(frameMat, opt.fs);
%		opt.polyOrder=-1;
%		opt.pitchRange=[40, 80];
%		opt.weight=0.9;
%		nshs=spec2nshs(powerSpec, opt, 1);

%	Roger Jang, 20101017

if nargin<1; selfdemo; return; end
[freqLen, frameNum]=size(specMat);
if freqLen==1, specMat=specMat(:); freqLen=frameNum; end	% Single frame
if nargin<2 | isempty(opt)
	opt.fs=8000;
	opt.polyOrder=6;			% Order for polynomial fitting for removing the trend of power spectrum
	opt.weight=0.8;			% Weighting factor for summation of power spectrum
	opt.pitchRange=[30, 84];
	opt.freq=opt.fs/(2*(freqLen-1))*(0:freqLen-1);
end
if nargin<3; plotOpt=0; end

specLen=length(opt.freq);
specMat2=specMat;
specTrend=specMat;
if opt.polyOrder>=0
	for i=1:frameNum
		specTrend(:,i)=polyval(polyfit((1:specLen)'/specLen, specMat(:,i), opt.polyOrder), (1:specLen)'/specLen);	% Divided by specLen to avoid bad scaling in polyfit
		specMat2(:,i)=specMat(:,i)-specTrend(:,i);		% Removing the trend
	end
end

freqRange=pitch2freq(opt.pitchRange);
%fprintf('freqRange=%s\n', mat2str(freqRange));
temp=find(opt.freq>=freqRange(1)); index1=temp(1);
temp=find(opt.freq<=freqRange(2)); index2=temp(end);

weightVec=(opt.weight).^(0:length(opt.freq)-1)';
nshs=nan*ones(length(opt.freq), frameNum);
for j=1:frameNum
	for i=index1:index2
		step=i-1;
		toBeAdded=specMat2(i:step:end, j);
		weight=weightVec(1:length(toBeAdded));
		nshs(i, j)=sum(toBeAdded.*weight)/sum(weight);
	end
end

if plotOpt 
	if frameNum==1
		plotNum=3;
		subplot(plotNum,1,1);
		plot(opt.freq, specMat, opt.freq, specTrend, 'r');
		set(gca, 'xlim', [-inf inf]);
		title('Power spectrum'); xlabel('Freq (Hz)');

		subplot(plotNum,1,2);
		plot(opt.freq, specMat2);
		set(gca, 'xlim', [-inf inf]);
		title('Power spectrum'); xlabel('Freq (Hz)');

		subplot(plotNum,1,3);
		plot(opt.freq, nshs); axis tight
		line(opt.freq(index1)*[1 1], [min(nshs), max(nshs)], 'color', 'r');
		line(opt.freq(index2)*[1 1], [min(nshs), max(nshs)], 'color', 'r');
		set(gca, 'xlim', [-inf inf]);
		title('Harmonic product spectrum'); xlabel('Freq (Hz)');
		[maxValue, maxIndex]=max(nshs);
		line(opt.freq(maxIndex), maxValue, 'color', 'r', 'marker', 'o');
	else
		subplot(1,2,1)
		imagesc(specMat); axis xy; axis image; colorbar
		xlabel('Frame index'); ylabel('Freq bin');
		title('Spectrum');
		subplot(1,2,2)
		imagesc(nshs); axis xy; axis image; colorbar
		xlabel('Frame index'); ylabel('Freq bin');
		title('NSHS map');
	end
end

% ====== Self demo
function selfdemo
waveFile='soo.wav';
startIndex=15000;
frameSize=512;
[y, opt.fs, nbits] = wavread(waveFile);
endIndex=startIndex+frameSize-1;
frame=y(startIndex:endIndex);
frame=frame.*hamming(frameSize);
opt.zeroPaddedFrameSize=16*frameSize;
opt.polyOrder=10;
opt.weight=0.8;
opt.pitchRange=[30, 84];
frame2=[frame; zeros(opt.zeroPaddedFrameSize-frameSize, 1)];		% zero padding
[mag, phase, opt.freq, pSpec]=fftOneSide(frame2, opt.fs, 0);
nshs=feval(mfilename, pSpec, opt, 1);
[maxValue, maxIndex]=max(nshs);
fprintf('Pitch freq = %f\n', opt.freq(maxIndex));
