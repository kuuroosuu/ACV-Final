function [pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt, plotOpt)
% pf2pitch: pf to pitch conversion
%	Usage: [pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt, plotOpt)

%	Roger Jang, 20070211

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

% ===== Find all local max
lmaxIndex1=find(localMax(pf));
lmaxValue1=pf(lmaxIndex1);

% ====== Within the right range
minIndex=floor(ptOpt.fs/pitch2freq(ptOpt.maxPitch));
maxIndex=ceil(ptOpt.fs/pitch2freq(ptOpt.minPitch));
index=find(minIndex<=lmaxIndex1 & lmaxIndex1<=maxIndex);
if isempty(index), pitch=0; clarity=0; pitchIndex=[]; return; end
lmaxValue2=lmaxValue1(index);
lmaxIndex2=lmaxIndex1(index);

% ====== Greater than ptOpt.clarityTh
index=find(lmaxValue2>=ptOpt.clarityTh);
lmaxIndex3=lmaxIndex2(index);
lmaxValue3=lmaxValue2(index);
if ~isempty(index)	% Find the left-most candidate
	% ====== Find the min. index of the local max
	[lmaxIndex4, index]=min(lmaxIndex3);
	lmaxValue4=lmaxValue3(index);
else			% Find the candidate with the max value
	[lmaxValue4, index]=max(lmaxValue2);
	lmaxIndex4=lmaxIndex2(index);
end

pitch=freq2pitch(ptOpt.fs/(lmaxIndex4-1));
clarity=lmaxValue4;
pitchIndex=lmaxIndex4;

if plotOpt
	plot(1:length(pf), pf);
	set(gca, 'xlim', [-inf inf]);
	axisLimit=axis;
	line(minIndex*[1 1], axisLimit(3:4), 'color', 'm');
	line(maxIndex*[1 1], axisLimit(3:4), 'color', 'm');
	line(minIndex*[1 1], [-1 1], 'color', 'm');
	line(maxIndex*[1 1], [-1 1], 'color', 'm');
	line(lmaxIndex1, lmaxValue1, 'color', 'g', 'linestyle', 'none', 'marker', 'o');
	line(lmaxIndex2, lmaxValue2, 'color', 'b', 'linestyle', 'none', 'marker', '^');
	line(lmaxIndex3, lmaxValue3, 'color', 'm', 'linestyle', 'none', 'marker', 'v');
	line(lmaxIndex4, lmaxValue4, 'color', 'r', 'linestyle', 'none', 'marker', 'square');
	line([1, length(pf)], ptOpt.clarityTh*[1, 1], 'color', 'r');
end

% ====== Self demo
function selfdemo
waveFile='greenOil.wav';
[y, fs, nbits]=wavReadInt(waveFile);
framedY=buffer2(y, 256, 0);
plotOpt=1;
ptOpt=ptOptSet(fs, nbits);
frame=framedY(:, 34);
pf=frame2nsdf(frame);
subplot(2,1,1); [pitch, clarity, ppcIndex]=feval(mfilename, pf, ptOpt, plotOpt);
frame=framedY(:, 250);
pf=frame2nsdf(frame);
subplot(2,1,2); [pitch, clarity, ppcIndex]=feval(mfilename, pf, ptOpt, plotOpt);