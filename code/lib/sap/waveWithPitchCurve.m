function waveWithPitchCurve(waveFile, outWaveFile, showPlot)
% waveWithPitchCurve: Read a wav file, perform pitch tracking, display the signals and the pitch in cooledit
%
%	Usage:
%		waveWithPitchCurve(waveFile, outWaveFile)
%
%	Example:
%		waveFile='whenWillYouBeBack.wav';
%		outWaveFile='whenWillYouBeBack_combined.wav';
%		waveWithPitchCurve(waveFile, outWaveFile, 1);

%	Category: Pitch tracking
%	Roger Jang, 20120706

if nargin<1, selfdemo; return; end
if nargin<2|isempty(outWaveFile);
	[parentDir, mainName, extName]=fileparts(waveFile);
	if isempty(parentDir)
		outWaveFile=[mainName, '_combined.wav'];
	else
		outWaveFile=[parentDir, '/', mainName, '_combined.wav'];
	end
end
if nargin<3, showPlot=0; end

waveFile='whenWillYouBeBack.wav';
wObj=myAudioRead(waveFile);
duration=length(wObj.signal)/wObj.fs;
wObj.signal=wObj.signal(1:duration*wObj.fs);
pfType=1;	% 0 for AMDF, 1 for ACF
ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
ptOpt.mainFun='maxPickingOverPf';
ptOpt.usePitchSmooth=0;
ptOpt.useVolThreshold=0;
ptOpt.useClarityThreshold=0;
[pitch, clarity]=pitchTracking(wObj, ptOpt, showPlot);

frameNum=length(pitch);
frameIndex=frame2sampleIndex(1:frameNum, ptOpt.frameSize, ptOpt.overlap);
newFrameIndex=[1:frameIndex(1)-1, frameIndex, frameIndex(end)+1:length(wObj.signal)];
newPitch=[zeros(1, frameIndex(1)-1), pitch, zeros(1, length(wObj.signal)-frameIndex(end))];

newPitch(newPitch==0)=nan;
minPitch=nanmin(newPitch);
newPitch(isnan(newPitch))=minPitch;
scaledPitch=0.8*(newPitch-min(newPitch))/(max(newPitch)-min(newPitch));	% The max pitch is 0.8 for display within cooledit
densePitch=interp1(newFrameIndex, scaledPitch, 1:length(wObj.signal));
wObj2=wObj;
wObj2.signal=[wObj.signal(:), densePitch(:)];	% Channel 1 is the original signals, channel 2 is the scaled pitch
fprintf('Saving %s...\n', outWaveFile);
wObj2file(wObj2, outWaveFile);
if showPlot
	dos(['start ', outWaveFile]);		% Start the wav file using cooledit
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);