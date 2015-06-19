function [pitch, clarity, pfMat, ptOpt, evalCount, maxPitchDiff]=pitchTrackingForcedSmooth(wObj, ptOpt, showPlot)
%pitchTrackingForcedSmooth: Pitch tracking which uses incremental weights to get a smooth pitch curve
%
%	Usage:
%		[pitch, clarity]=pitchTrackingForcedSmooth(waveObj, ptOpt, showPlot);
%
%	Example:
%		waveFile='test01.wav';
%		aObj=myAudioRead(waveFile);
%		ptOpt=ptOptSet(aObj.fs, aObj.nbits);
%		ptOpt.pitchDiffTh=6;
%		[pitch, clarity, pfMat, ptOpt, evalCount, maxPitchDiff]=pitchTrackingForcedSmooth(aObj, ptOpt, 1);

if nargin<1, selfdemo; return; end
if nargin<2, showPlot=1; end

if nargin<1, selfdemo; return; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name
if nargin<2 | isempty(ptOpt), ptOpt=ptOptSet(wObj.fs, wObj.nbits); end
if nargin<3, showPlot=0; end
if wObj.fs~=ptOpt.fs | wObj.nbits~=ptOpt.nbits
	fprintf('wObj.fs=%d, ptOpt.fs=%d, wObj.nbits=%d, ptOpt.nbits=%d\n', wObj.fs, ptOpt.fs, wObj.nbits, ptOpt.nbits);
	error('ptOpt is not conformed to the given wave object!');
end

% ====== Define parameters
pitchDiffTh=ptOpt.pitchDiffTh;
roundNum=ptOpt.roundNum;
searchMethod=ptOpt.searchMethod;	% 'binarySearch' or 'linearSearch'
leftBound=ptOpt.leftBound; rightBound=ptOpt.rightBound;	% for binary search only
boundRangeTh=ptOpt.boundRangeTh;

maxPitchDiff=inf;
evalCount=0;
switch(searchMethod)
	case 'linearSearch'
		while maxPitchDiff>pitchDiffTh
			ptOpt.indexDiffWeight=ptOpt.indexDiffWeight+100;
			[pitch, clarity, pfMat, maxPitchDiff]=pitchTracking(wObj, ptOpt); evalCount=evalCount+1;
			pitchData(evalCount).indexDiffWeight=ptOpt.indexDiffWeight;
			pitchData(evalCount).pitch=pitch;
			legendStr{evalCount}=num2str(ptOpt.indexDiffWeight);
			fprintf('\tindexDiffWeight=%f, evalCount=%d, maxPitchDiff=%f\n', ptOpt.indexDiffWeight, evalCount, maxPitchDiff);
		end
	case 'binarySearch'
		% ====== Find the right rightBound which make the pitch diff less than a threshold
		fprintf('Initial bounds: [leftBound, rightBound]=[%g, %g], pitchDiffTh=%g\n', leftBound, rightBound, pitchDiffTh);
		ptOpt.indexDiffWeight=leftBound; [pitchLeft, clarity, pfMat, maxPitchDiff]=pitchTracking(wObj, ptOpt); evalCount=evalCount+1;
		while maxPitchDiff>pitchDiffTh
			ptOpt.indexDiffWeight=rightBound;
			[pitchRight, clarity, pfMat, maxPitchDiff]=pitchTracking(wObj, ptOpt); evalCount=evalCount+1;
			fprintf('Finding the right rightBound: leftBound=%g, rightBound=%g, maxPitchDiff=%g\n', leftBound, rightBound, maxPitchDiff);
			if maxPitchDiff<pitchDiffTh, break; end
			leftBound=rightBound; pitchLeft=pitchRight;
			rightBound=rightBound*2;
		end
		% ====== Shorten the bound by binary search
		for i=1:roundNum
			center=(leftBound+rightBound)/2;
			ptOpt.indexDiffWeight=center;
			[pitchCenter, clarity, pfMat, maxPitchDiff]=pitchTracking(wObj, ptOpt); evalCount=evalCount+1;
			if maxPitchDiff>pitchDiffTh
				leftBound=center; pitchLeft=pitchCenter;
			else
				rightBound=center; pitchRight=pitchCenter;
			end
			fprintf('Tighten the bound: i=%d, leftBound=%g, rightBound=%g, boundRange=%g, maxPitchDiff=%g at %g\n', i, leftBound, rightBound, rightBound-leftBound, maxPitchDiff, center);
			if rightBound-leftBound<boundRangeTh, break; end
		end
		ptOpt.indexDiffWeight=rightBound; pitch=pitchRight;
	otherwise
		error('Unknown method!');
end

if showPlot
	[pitch, clarity, pfMat]=pitchTracking(wObj, ptOpt, 1);	% Deliver the final plot
	if strcmp(searchMethod, 'linearSearch')
		curveNum=evalCount-1;
		allPitch=cat(1, pitchData.pitch)';
		allPitch(allPitch==0)=nan;
		for i=1:curveNum, allPitch(:,i)=allPitch(:,i)+0.3*(i-1); end
		figure; plot(allPitch);
		legend(legendStr, 'location', 'best');
		title('Pitch curves vs indexDiffWeight (for AMDF)');
	end
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
