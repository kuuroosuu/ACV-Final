function [pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt, plotOpt)
% pf2pitch: pf to pitch conversion
%	Usage: [pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt, plotOpt)

%	Roger Jang, 20070211

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

switch(ptOpt.pfType)
	case 0	% AMDptOpt
		optFun='min';
		localOptFun='localMin';
	case 1	% ACF
		optFun='max';
		localOptFun='localMax';
	otherwise
		error('Unknown pfType!');
end

% ===== Find all local optima
localOptIndex1=find(feval(localOptFun, pf));
localOptValue1=pf(localOptIndex1);

% ====== Within the right range
minIndex=floor(ptOpt.fs/pitch2freq(ptOpt.maxPitch));
maxIndex=ceil(ptOpt.fs/pitch2freq(ptOpt.minPitch));
index=find(minIndex<=localOptIndex1 & localOptIndex1<=maxIndex);
if isempty(index), pitch=0; clarity=0; pitchIndex=[]; return; end
localOptValue2=localOptValue1(index);
localOptIndex2=localOptIndex1(index);

% ====== Find the optimum PF value
[optPfValue, localIndex]=feval(optFun, localOptValue2);
pitchIndex=localOptIndex2(localIndex);
pitch=freq2pitch(ptOpt.fs/(pitchIndex-1));
switch(ptOpt.pfType)
	case 0	% AMDF
		clarity=1-optPfValue/max(pf(minIndex:maxIndex));		
	case 1	% ACF
		clarity=optPfValue/pf(1);
	otherwise
		error('Unknown pfType!');
end

if plotOpt
	plot(1:length(pf), pf, '.-');
	set(gca, 'xlim', [-inf inf]);
	axisLimit=axis;
	line(minIndex*[1 1], axisLimit(3:4), 'color', 'm');
	line(maxIndex*[1 1], axisLimit(3:4), 'color', 'k');
	line(localOptIndex1, localOptValue1, 'color', 'g', 'linestyle', 'none', 'marker', 'o');
	line(localOptIndex2, localOptValue2, 'color', 'b', 'linestyle', 'none', 'marker', '^');
	line(pitchIndex, optPfValue, 'color', 'r', 'linestyle', 'none', 'marker', 'square');
	legend('AMDF', 'Lower bound', 'Upper bound', 'All local min', 'Within-bound local min', 'Selected');
end

% ====== Self demo
function selfdemo
waveFile='greenOil.wav';
[y, fs, nbits]=wavReadInt(waveFile);
framedY=buffer2(y, 256, 0);
plotOpt=1;
ptOpt=ptOptSet(fs, nbits, 0);
frame=framedY(:, 34);
pf=frame2amdf(frame, length(frame), 2);
subplot(2,1,1); [pitch, clarity, ppcIndex]=feval(mfilename, pf, ptOpt, plotOpt);
frame=framedY(:, 250);
pf=frame2amdf(frame, length(frame), 2);
subplot(2,1,2); [pitch, clarity, ppcIndex]=feval(mfilename, pf, ptOpt, plotOpt);