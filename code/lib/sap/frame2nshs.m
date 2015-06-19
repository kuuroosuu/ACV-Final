function [nshs, freq]=frame2nshs(frameMat, opt, plotOpt)
% frame2nshs: Conversion from a frame (or frame matrix) to normalized sub-harmonic sum
%
%	Usage:
%		[nshs, freq]=frame2nshs(frameMat, opt, plotOpt)
%
%	For example:
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
%		[nshs, freq]=frame2nshs(frameMat, [], 1);

%	Roger Jang, 20101017

if nargin<1; selfdemo; return; end
[frameSize, frameNum]=size(frameMat);
if frameSize==1, frameMat=frameMat(:); frameSize=frameNum; end	% Single frame
if nargin<2 | isempty(opt)
	opt.fs=8000;
	opt.zeroPaddedFrameSize=16*length(frameSize);
	opt.polyOrder=10;			% Order for polynomial fitting for removing the trend of power spectrum
	opt.weight=0.8;				% Weighting factor for summation of power spectrum
	opt.pitchRange=[30, 84];
end
if nargin<3; plotOpt=0; end

frameMat2=frameMat;
for i=1:frameNum
	frameMat2(:,i)=frameMat(:,i).*hamming(frameSize);	% Hamming windowing
end

if opt.zeroPaddedFrameSize>frameSize
	frameMat2=[frameMat2; zeros(opt.zeroPaddedFrameSize-frameSize, frameNum)];		% zero padding
end

[mag, phase, freq, specMat]=fftOneSide(frameMat2, opt.fs, 0);
opt.freq=freq;
[nshs, specMat2, specTrend]=spec2nshs(specMat, opt);

if plotOpt 
	if frameNum==1
		frame=frameMat(:);
		plotNum=4;
		subplot(plotNum,1,1);
		plot(1:length(frame), frame, '.-'); axis tight;
		title('Frame'); xlabel('Samples');

		subplot(plotNum,1,2);
		plot(freq, specMat, freq, specTrend, 'r');
		set(gca, 'xlim', [-inf inf]);
		title('Power spectrum'); xlabel('Freq (Hz)');

		subplot(plotNum,1,3);
		plot(freq, specMat2);
		set(gca, 'xlim', [-inf inf]);
		title('Power spectrum'); xlabel('Freq (Hz)');

		subplot(plotNum,1,4);
		freqRange=pitch2freq(opt.pitchRange);
		temp=find(opt.freq>=freqRange(1)); index1=temp(1);	% For plotting
		temp=find(opt.freq<=freqRange(2)); index2=temp(end);	% For plotting
		plot(freq, nshs); axis tight
		line(freq(index1)*[1 1], [min(nshs), max(nshs)], 'color', 'r');
		line(freq(index2)*[1 1], [min(nshs), max(nshs)], 'color', 'r');
		set(gca, 'xlim', [-inf inf]);
		title('Harmonic product spectrum'); xlabel('Freq (Hz)');
		[maxValue, maxIndex]=max(nshs);
		line(freq(maxIndex), maxValue, 'color', 'r', 'marker', 'o');
	else
		subplot(1,3,1)
		imagesc(specMat); axis xy; axis image; colorbar
		title('Original power spectrum');
		subplot(1,3,2)
		imagesc(specMat2); axis xy; axis image; colorbar
		title('Formant-removed power spectrum');
		subplot(1,3,3)
		imagesc(nshs); axis xy; axis image; colorbar
		title('NSHS map');
		specMat(specMat<-10000)=nan;
		figure; mesh(specMat);
		figure; mesh(specMat2);
	end
end

% ====== Self demo
function selfdemo
waveFile='\dataSet\dean\201005-mirRecording\Record_Wave/陳亮宇_1/中華民國國歌_不詳_0_s.wav';
startIndex=25098;
frameSize=512;
waveFile='soo.wav';
startIndex=15000;
frameSize=512;
[y, opt.fs, nbits] = wavread(waveFile);
endIndex=startIndex+frameSize-1;
frame=y(startIndex:endIndex);
opt.zeroPaddedFrameSize=16*frameSize;
opt.polyOrder=10;
opt.weight=0.8;
opt.pitchRange=[30, 84];
[nshs, freq]=feval(mfilename, frame, opt, 1);
[maxValue, maxIndex]=max(nshs);
fprintf('Pitch freq = %f\n', freq(maxIndex));
