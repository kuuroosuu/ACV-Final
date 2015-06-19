function [pitch, clarity, pfMat, maxPitchDiff]=pitchTrack(wObj, ptOpt, showPlot)
% pitchTrack: pitch tracking for a given file
%
%	Usage:
%		[pitch, clarity]=pitchTrack(waveFile, ptOpt, showPlot)
%
%	Example:
%		% === Example of singing input
%		waveFile='yankee_doodle.wav';
%		waveFile='twinkle_twinkle_little_star.wav';
%	%	waveFile='10LittleIndians.wav';	% Noisy!
%		wObj=waveFile2obj(waveFile);
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits);
%		pitch=pitchTrack(wObj, ptOpt, 1);
%		% === Example of speech
%		waveFile='what_movies_have_you_seen_recently.wav';
%		waveFile='但使龍城飛將在.wav';
%		wObj=waveFile2obj(waveFile);
%		pfType=1;	% 0 for AMDF, 1 for ACF
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
%		ptOpt.frameSize=640; ptOpt.overlap=640-160;	% frameSize 是 HTK 的兩倍，frame rate=100, 搭配語音辨識使用
%		figure; [pitch, clarity, pfMat]=pitchTrack(wObj, ptOpt, 1);
%		% === Example with mainFun = 'maxPickingOverPf'
%		waveFile='twinkle_twinkle_little_star.wav';
%		wObj=waveFile2obj(waveFile);
%		pfType=1;	% 0 for AMDF, 1 for ACF
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
%		ptOpt.mainFun='maxPickingOverPf';
%		showPlot=1;
%		figure; [pitch, clarity]=pitchTrack(wObj, ptOpt, showPlot);
%		% === Another example without pitch groundtruth
%		waveFile='soo.wav';
%		wObj=waveFile2obj(waveFile);
%		pfType=1;		% 0 for AMDF, 1 for ACF
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
%		ptOpt.mainFun='maxPickingOverPf';
%		ptOpt.usePitchSmooth=0;
%		showPlot=1;
%		figure; [pitch, clarity]=pitchTrack(wObj, ptOpt, showPlot);

%	Roger Jang, 20090810, 20121008

if nargin<1, selfdemo; return; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name
if nargin<2 | isempty(ptOpt), ptOpt=ptOptSet(wObj.fs, wObj.nbits); end
if nargin<3, showPlot=0; end
if wObj.fs~=ptOpt.fs | wObj.nbits~=ptOpt.nbits
	error('wObj.fs=%d, ptOpt.fs=%d, , wObj.nbits=%d, ptOpt.nbits=%d, ptOpt is not conformed to the given wave object!', wObj.fs, ptOpt.fs, wObj.nbits, ptOpt.nbits);
end

% Make wObj.signal integer for ptByDpOverPfMex
if wObj.amplitudeNormalized
	wObj.signal=wObj.signal*2^(wObj.nbits-1);
	wObj.amplitudeNormalized=0;
end

if ptOpt.useWaveEnhancement
	fprintf('\t\t\tPerform speech enhancement...\n');
	wObj.signal=specSubtractMex(wObj.signal, wObj.fs);
end

if ptOpt.useHighPassFilter
	fprintf('\t\t\tPerform high-pass filtering...\n');
	cutOffFreq=100;		% Cutoff frequency
	filterOrder=5;		% Order of filter
	[b, a]=butter(filterOrder, cutOffFreq/(wObj.fs/2), 'high');
	wObj.signal=filter(b, a, wObj.signal);
end

if ptOpt.useLowPassFilter
	fprintf('\t\t\tPerform low-pass filtering...\n');
	cutOffFreq=500;		% Cutoff frequency
	filterOrder=5;		% Order of filter
	[b, a]=butter(filterOrder, cutOffFreq/(wObj.fs/2), 'low');
	wObj.signal=filter(b, a, wObj.signal);
end

y=wObj.signal; fs=wObj.fs; nbits=wObj.nbits;
% ====== Pitch tracking
switch(ptOpt.mainFun)
	case 'dpOverPfMat'
	%	[pitch, clarity, pfMat]=ptByDpOverPfMex(wObj, ptOpt);
		[pitch, clarity, pfMat, dpPath]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
	case 'maxPickingOverPf'
		frameMat=enframe(y, ptOpt.frameSize, ptOpt.overlap);
		frameMat=frameZeroMean(frameMat, ptOpt.zeroMeanPolyOrder);
		frameNum=size(frameMat, 2);
		pitch=zeros(1, frameNum);
		clarity=zeros(1, frameNum);
		for i=1:frameNum
			frame=frameMat(:, i);
			[pitch(i), clarity(i), pitchIndex, pf]=frame2pitch2(frame, ptOpt);	% This should by changed to frame2pitch()!
			lMaxCount(i)=sum(localMax(pf));
		end
		pitch2=pitch;	% Keep pitch2 for plotting
		% ===== Remove pitch with low volume
		if ptOpt.useVolThreshold
			epdPrm=endPointDetect('defaultOpt');
			epdPrm.frameDuration=ptOpt.frameSize/ptOpt.fs;
			epdPrm.overlapDuration=ptOpt.overlap/ptOpt.fs;	
			epdPrm.method='vol';
			[epInSampleIndex, epInFrameIndex, segment, zeroOneVec, others]=endPointDetect(wObj, epdPrm);
			volume=others.volume; volTh=others.volTh;
%			[epInSampleIndex, epInFrameIndex, segment, zeroOneVec, volume, volTh]=epdByVol(y, ptOpt.fs, ptOpt.nbits, epdPrm);
			pitch=pitch.*zeroOneVec;
		end
		% ===== Remove pitch with low clarity
		if ptOpt.useClarityThreshold
			clarityTh=ptOpt.clarityRatio*max(clarity);
			pitch(clarity<clarityTh)=0;
		end
		% ====== Remove pitch with high pf local max count
		if ptOpt.localMaxCountThresholding
			lMaxCountTh=ptOpt.localMaxCountRatio*max(lMaxCount);
			pitch(lMaxCount>lMaxCountTh)=0;
		end
		% ====== Smooth pitch
		if ptOpt.usePitchSmooth
			pitch=pitchSmooth(pitch);
		end
	otherwise
		error(sprintf('Unknown method "%s" in pitchTrack()!', ptOpt.mainFun));
end

% ====== SU/V detection via HMM
if ptOpt.useHmm4suvDetection
	fprintf('\t\t\tPerform SU/V detection via HMM...\n');
	load pitchExist/gmmData.mat
	load pitchExist/transProb.mat
	transLogProb=log(transProb);
	gmmSet=gmmData(index);
	suvParam=suvParamSet;
	suv=wave2suv(wObj, suvParam, gmmSet, transLogProb);
	pitch=pitch.*(suv'-1);
end

if nargout>=4
	segment=segmentFind(pitch);
	segmentNum=length(segment);
	segmentPitchDiff=zeros(1,segmentNum);
	for i=1:segmentNum
		if segment(i).duration>1
			segmentPitchDiff(i)=max(abs(diff(pitch(segment(i).begin:segment(i).end))));
		end
	end
	maxPitchDiff=max(segmentPitchDiff);
end

%{
frameMat=enframe(y, ptOpt.frameSize, ptOpt.overlap);
frameMat=frameZeroMean(frameMat, 2);
frameNum=size(frameMat, 2);
clarity=frame2clarity(frameMat, fs, 'acf', 1);
%}
[maxClarity, maxIndex]=max(clarity);
clarityTh=ptOpt.clarityRatio*maxClarity;
%fprintf('maxClarity=%f, clarityTh=%f\n', max(clarity), clarityTh);
%pitch=pitch.*(clarity>clarityTh);

if showPlot
	plotNum=2;
	if strcmp(ptOpt.mainFun, 'dpOverPfMat'), plotNum=plotNum+1; end
	if ptOpt.useVolThreshold, plotNum=plotNum+1; end
	if ptOpt.useClarityThreshold, plotNum=plotNum+1; end
	if ptOpt.localMaxCountThresholding, plotNum=plotNum+1; end
	plotId=0;
	frameNum=floor((length(y)-ptOpt.overlap)/(ptOpt.frameSize-ptOpt.overlap));
	frameTime=frame2sampleIndex(1:frameNum, ptOpt.frameSize, ptOpt.overlap)/fs;

	plotId=plotId+1; subplot(plotNum, 1, plotId);
	time=(1:length(y))/fs;
	plot(time, y); axis([min(time) max(time) 2^nbits/2*[-1 1]]);
	if ~isempty(wObj.file) title(strPurify4label(sprintf('Waveform of %s', wObj.file))); end

	if strcmp(ptOpt.mainFun, 'dpOverPfMat')
		plotId=plotId+1; subplot(plotNum, 1, plotId);
		pfPlot(pfMat, frameTime, dpPath);
		title('PF matrix (white dots: DP path, black dots: Pitch after all kinds of thresholding/smoothing');
		pp=round(fs./pitch2freq(pitch));
		for i=1:frameNum, line(frameTime(i), pp(i), 'color', 'k', 'marker', '.'); end
	end

	plotId=plotId+1; subplot(plotNum, 1, plotId);
	cPitch=pitch; cPitch(cPitch==0)=nan;
	tPitch=nan*cPitch;
	pvFile=[wObj.file(1:end-3), 'pv'];
	if exist(pvFile), wObj.tPitch=asciiRead(pvFile); end
	if isfield(wObj, 'tPitch')
		tPitch=wObj.tPitch;
		tPitch(tPitch==0)=nan;
		plot(frameTime, tPitch, 'r-o', frameTime, cPitch, 'k.-');	% Do not change the order!
		title('Desired (circle) and computed (dot) pitch');
	%	legend('Computed pitch', 'Desired pitch', 'location', 'northOutside', 'orientation', 'horizontal');
	else
		plot(frameTime, cPitch, 'k.-');
		title('Computed pitch');
	%	legend('Computed pitch', 'location', 'northOutside', 'orientation', 'horizontal');
	end
%	line([min(frameTime), max(frameTime)], ptOpt.minPitch*[1 1], 'color', 'm');
%	line([min(frameTime), max(frameTime)], ptOpt.maxPitch*[1 1], 'color', 'm');
%	axis([min(frameTime) max(frameTime) ptOpt.minPitch-1 ptOpt.maxPitch+1]);
%	axis([min(frameTime) max(frameTime) -inf inf]);
	set(gca, 'xlim', [min(frameTime) max(frameTime)]); axis tight
	ylabel('Pitch (semitone)'); grid on

	if ptOpt.useVolThreshold
		plotId=plotId+1; subplot(plotNum, 1, plotId);
		frameTime=frame2sampleIndex(1:length(pitch), ptOpt.frameSize, ptOpt.overlap)/ptOpt.fs;
		% EPD, which should have been done in ptByDpOverPfMex!!!
		epdPrm=endPointDetect('defaultOpt');
		epdPrm.frameDuration=ptOpt.frameSize/ptOpt.fs;
		epdPrm.overlapDuration=ptOpt.overlap/ptOpt.fs;	
		epdPrm.method='vol';
		[epInSampleIndex, epInFrameIndex, segment, zeroOneVec, others]=endPointDetect(wObj, epdPrm);		
		volume=others.volume; volTh=others.volTh;
%		[epInSampleIndex, epInFrameIndex, segment, zeroOneVec, volume, volTh]=epdByVol(y, ptOpt.fs, ptOpt.nbits, epdPrm);
		pitch=pitch.*zeroOneVec;
		plot(frameTime, volume, '.-b', [min(frameTime), max(frameTime)], volTh*[1 1], 'r'); title('Volume');
	%	line([min(frameTime), max(frameTime)], volTh*[1 1], 'color', 'r');
		for i=1:length(segment)
			line(frameTime(segment(i).beginFrame)*[1 1], [0, max(volume)], 'color', 'g');
			line(frameTime(segment(i).endFrame)*[1 1], [0, max(volume)], 'color', 'm');
		end
		set(gca, 'xlim', [min(frameTime) max(frameTime)]);
	end

	if ptOpt.useClarityThreshold
		plotId=plotId+1; subplot(plotNum, 1, plotId);
		plot(frameTime, clarity, '.-');
		line([min(frameTime), max(frameTime)], max(clarity)*ptOpt.clarityRatio*[1 1], 'color', 'r');
		set(gca, 'xlim', [min(frameTime) max(frameTime)]);
		xlabel('Time (sec)'); title('Clarity'); grid on
	end
	% ====== Buttons for playback
	pitchObj.signal=pitch; pitchObj.frameRate=fs/(ptOpt.frameSize-ptOpt.overlap);
	if isfield(wObj, 'tPitch')
		pitchObj2.signal=tPitch; pitchObj2.frameRate=fs/(ptOpt.frameSize-ptOpt.overlap);
		buttonH=wavePitchPlayButton(wObj, pitchObj, pitchObj2);
		set(buttonH(end), 'string', 'Play GT pitch')
	else
		buttonH=wavePitchPlayButton(wObj, pitchObj);
	end
%	set(gcf, 'name', sprintf('PT using ptByDpOverPfMex, with pfWeight=%d and indexDiffWeight=%d', ptOpt.pfWeight, ptOpt.indexDiffWeight));
	set(gcf, 'name', sprintf('PT using mainFun=%s', ptOpt.mainFun));
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
