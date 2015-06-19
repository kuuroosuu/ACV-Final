function onset=odByVol(wObj, odPrm, plotOpt)
% tpByVol: Onset detection of tapping by volume
%	Usage: onset=odByVol(y, fs, odPrm, plotOpt)
%
%	For example:
%		waveFile='tapping.wav';
%		[y, fs, nbits, opts, cueLabel]=wavReadInt(waveFile);
%		plotOpt=1;
%		odPrm=odPrmSet;
%		onset=odByVol(y, fs, odPrm, plotOpt);
%		subplot(2,1,1);
%		axisLimit=axis;
%		% Display the human-transcribed cue labels
%		line(cueLabel/fs, axisLimit(4)*ones(length(cueLabel),1), 'color', 'r', 'marker', 'v', 'linestyle', 'none'); 

%	Roger Jang, 20090607

if nargin<1, selfdemo; return; end
if nargin<3, odPrm=odPrmSet; end
if nargin<4, plotOpt=0; end

y=wObj.signal; fs=wObj.fs;

opt=wave2volume('defaultOpt');
opt.frameSize=round(odPrm.frameSize*fs);
opt.overlap=round((odPrm.frameSize-odPrm.frameShift)*fs);
opt.frame2volumeOpt.method='absSum';

if odPrm.useHighPassFilter
	[b, a]=butter(odPrm.filterOrder, odPrm.cutOffFreq/(fs/2), 'high');
	y=filter(b, a, y);
end

% ====== Step 1: Apply a volume threshold to have onset candidates
volume=wave2volume(wObj, opt);
sortedVol=sort(volume);
vol1=sortedVol(round(length(volume)*4/5));
vol2=sortedVol(round(length(volume)*99/100));
volTh=(vol2-vol1)*odPrm.volRatio+vol1;
frameNum=length(volume);
onset1=volume>volTh & localMax(volume);

% ====== Step 2: Apply a moving window to select the right onset
onset2=onset1;
halfWinWidth=fs/odPrm.maxTappingPerSec/(opt.frameSize-opt.overlap);
index=find(onset2);
for i=index
	startIndex=max(i-round(halfWinWidth), 1);
	endIndex=min(i+round(halfWinWidth), frameNum);
	[junk, maxIndex]=max(volume(startIndex:endIndex));
	if maxIndex+startIndex-1~=i
		onset2(i)=0;
	end
end
onset=frame2sampleIndex(find(onset2), opt.frameSize, opt.overlap);

% Plotting
if plotOpt
	sampleTime=(1:length(y))/fs;
	frameTime=frame2sampleIndex(1:length(volume), opt.frameSize, opt.overlap)/fs;
	subplot(2,1,1);
	plot(sampleTime, y); set(gca, 'xlim', [-inf inf]);
	% Display the detected tapping
	axisLimit=axis;
	line(onset/fs, axisLimit(3)*ones(length(onset),1), 'color', 'k', 'marker', '^', 'linestyle', 'none');
	xlabel('Time (sec)'); ylabel('Waveform');
	subplot(2,1,2);
	plot(frameTime, volume, '.-'); set(gca, 'xlim', [-inf inf]);
	xlabel('Time (sec)'); ylabel('Volume');
	line([frameTime(1), frameTime(end)], volTh*[1 1], 'color', 'r');
	line([frameTime(1), frameTime(end)], vol1*[1 1], 'color', 'r');
	line([frameTime(1), frameTime(end)], vol2*[1 1], 'color', 'r');
	line(frameTime(onset1), volume(onset1), 'marker', '.', 'color', 'g', 'linestyle', 'none');
	line(frameTime(onset2), volume(onset2), 'marker', '^', 'color', 'k', 'linestyle', 'none');
end

% ====== Self demo
function selfdemo
waveFile='tapping.wav';
[y, fs, nbits, opts, cueLabel]=wavReadInt(waveFile);
wObj=waveFile2obj(waveFile);
plotOpt=1;
odPrm=odPrmSet;
onset=feval(mfilename, wObj, odPrm, plotOpt);
subplot(2,1,1);
axisLimit=axis;
% Display the human-transcribed cue labels
line(cueLabel/fs, axisLimit(4)*ones(length(cueLabel),1), 'color', 'r', 'marker', 'v', 'linestyle', 'none'); 
