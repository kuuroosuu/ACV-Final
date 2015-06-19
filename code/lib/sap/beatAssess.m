function [score, avgDist]=beatAssess(wObj, opt, showPlot)
%beatAssess: Scoring of beat synchronization
%
%	Example:
%		waveFile='D:\dataSet\metronome\set03_iPod\13_50.wav';
%		wObj=waveFile2obj(waveFile);
%		opt=beatAssess('defaultOpt');
%		showPlot=1;
%		[score, avgDist]=beatAssess(wObj, opt, showPlot);

%	Roger Jang, 20130430

if nargin<1, selfdemo; return; end
if nargin==1 && ischar(wObj) && strcmpi(wObj, 'defaultOpt')
	score.numLeadingTick=4;
	score.pureTickDuration=3;	% sec
	score.cutRatio4tickSpectrum=0.05;
	score.btOpt=btOptSet;
	return
end
if nargin<3, showPlot=0; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name

%% Identify cutIndex to separate the first few clean ticks
if showPlot, figure; end
[ncObj, rawNc, smoothedNc, localMeanNc, msSpec]=wave2osc(wObj, opt.btOpt.oscOpt, showPlot);
localMaxIndex=find(localMax(ncObj.signal));
for i=1:length(localMaxIndex)
	index=localMaxIndex(i);
	if showPlot
		line(ncObj.time(index), ncObj.signal(index), 'color', 'r', 'marker', '.');
	end
end
if isfield(opt, 'pureTickDuration')	% Use this to define cutIndex directly
	cutIndex=opt.pureTickDuration*wObj.fs;
else	% Identify the cutIndex
	maxPdf=max(ncObj.signal);
	aboveAvgIndex=find(ncObj.signal>maxPdf/2);
	candidateIndex=intersect(localMaxIndex, aboveAvgIndex);
	cutIndex=round(mean(ncObj.time(candidateIndex(opt.numLeadingTick:opt.numLeadingTick+1))*wObj.fs));
end
if showPlot
	subplot(414); axisLimit=axis;
	line(cutIndex/wObj.fs*[1 1], axisLimit(3:4), 'color', 'm', 'linewidth', 3);
end
%% Use signals before the cutIndex to identify beatPeriod and unwanted frequency bins
wObj2=wObj; wObj2.signal=wObj2.signal(1:cutIndex);	% wObj2 is pure ticks
if showPlot, figure; plot((1:cutIndex)/wObj2.fs, wObj2.signal); end
% ====== Plot novelty curve
if showPlot, figure; end
[ncObj, rawNc, smoothedNc, localMeanNc, msSpec]=wave2osc(wObj2, opt.btOpt.oscOpt, showPlot);
% ====== Perform beat tracking to find the IBI
opt.btOpt.type='constant';
cBeat=beatTracking(wObj2, opt.btOpt);
ibi=diff(cBeat);
tickPeriod=mean(ibi(2:end));	% Assuming the first tick is not stable
% ====== Find the frequency bins not to be used in computing NC
sumSpec=sum(msSpec,2);
th=(max(sumSpec)-min(sumSpec))*0.05+min(sumSpec);
deleteIndex=find(sumSpec>th);
if showPlot
	figure; plot(sumSpec, '.-');
	axisLimit=axis;
	line(axisLimit(1:2), th*[1 1], 'color', 'r');
end
%% ====== Compute NC for the whole recording
opt.btOpt.deleteIndex=deleteIndex;
if showPlot, figure; end
[ncObj, rawNc, smoothedNc, localMeanNc, msSpec]=wave2osc(wObj, opt.btOpt.oscOpt, showPlot);
[~, mainName]=fileparts(wObj.file);
if showPlot,
	subplot(411), title(sprintf('Mel-band spectrogram: %s.wav', strPurify4label(mainName)));
end
% ====== Find tick positions
nextBeatPos=cBeat(end)+tickPeriod;
i=1; tickTime(i)=nextBeatPos;
while nextBeatPos<=(length(wObj.signal)-1)/wObj.fs
	nextBeatPos=nextBeatPos+tickPeriod;
	i=i+1; tickTime(i)=nextBeatPos;
end
% ====== Plot the tick position
if showPlot
	subplot(414);
	axisLimit=axis;
	for i=1:length(cBeat)
		line(cBeat(i)*[1 1], axisLimit(3:4), 'color', 'g');
	end
	line(cutIndex/wObj.fs*[1 1], axisLimit(3:4), 'color', 'm', 'linewidth', 3);
	for i=1:length(tickTime)
		line(tickTime(i)*[1 1], axisLimit(3:4), 'color', 'r');
	end
end
%% ====== Compute distance and score
localMaxIndex=find(localMax(ncObj.signal));
localMaxTime=ncObj.time(localMaxIndex);
dist=distPairwise(localMaxTime, tickTime);
[minDistVec, nearestOnsetIndex]=min(dist);
avgDist=mean(minDistVec);
score=100/(1+10*avgDist);

if showPlot
	for i=1:length(nearestOnsetIndex)
		index=localMaxIndex(nearestOnsetIndex(i));
		line(ncObj.time(index), ncObj.signal(index), 'color', 'r', 'marker', '.'); 
	end
end
% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);