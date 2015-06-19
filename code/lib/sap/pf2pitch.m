function [pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt, plotOpt)
% pf2pitch: PF (pitch function) to pitch conversion
%
%	Usage:
%		[pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt, plotOpt)
%
%	Example:
%		% === pfType=0 (AMDF)
%		waveFile='greenOil.wav';
%		[y, fs, nbits]=wavReadInt(waveFile);
%		frameMat=buffer2(y, 256, 0);
%		pfType=0;
%		ptOpt=ptOptSet(fs, nbits, pfType);
%		frameId=34;
%		frame=frameMat(:, frameId); pf=frame2pf(frame, pfType);
%		subplot(2,1,1); [pitch, clarity, ppcIndex]=pf2pitch(pf, ptOpt, 1); grid on
%		title(sprintf('pfType=%d, frameId=%d\n', pfType, frameId));
%		frameId=250;
%		frame=frameMat(:, frameId); pf=frame2pf(frame, pfType);
%		subplot(2,1,2); [pitch, clarity, ppcIndex]=pf2pitch(pf, ptOpt, 1); grid on
%		title(sprintf('pfType=%d, frameId=%d\n', pfType, frameId));
%		% === pfType=1 (ACF)
%		pfType=1;
%		ptOpt=ptOptSet(fs, nbits, pfType);
%		frameId=34;
%		frame=frameMat(:, frameId); pf=frame2pf(frame, pfType);
%		figure;
%		subplot(2,1,1); [pitch, clarity, ppcIndex]=pf2pitch(pf, ptOpt, 1); grid on
%		title(sprintf('pfType=%d, frameId=%d\n', pfType, frameId));
%		frameId=250;
%		frame=frameMat(:, frameId); pf=frame2pf(frame, pfType);
%		subplot(2,1,2); [pitch, clarity, ppcIndex]=pf2pitch(pf, ptOpt, 1); grid on
%		title(sprintf('pfType=%d, frameId=%d\n', pfType, frameId));
%
%	Roger Jang, 20070211, 20100617

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

switch(ptOpt.pfType)
	case 0	% AMDF
		optFun='min';
		localOptFun='localMin';
	case 1	% ACF
		optFun='max';
		localOptFun='localMax';
	case 2	% NSDF
		optFun='max';
		localOptFun='localMax';
	case 3	% HPS
		optFun='max';
		localOptFun='localMax';
	otherwise
		error(sprintf('Unknown pfType=%d!', ptOpt.pfType));
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
if ptOpt.useParabolicFit
	deviation=optimViaParabolicFit(pf(pitchIndex-1:pitchIndex+1));
	finetunedPitchIndex=pitchIndex+deviation;
	pitch=freq2pitch(ptOpt.fs/(finetunedPitchIndex-1));
else
	pitch=freq2pitch(ptOpt.fs/(pitchIndex-1));
end

switch(ptOpt.pfType)
	case 0	% AMDF
		clarity=1-optPfValue/max(pf(minIndex:maxIndex));		
	case 1	% ACF
		clarity=optPfValue/pf(1);
	case 2	% NSDF
		clarity=pf(1);
	case 3	% HPS
		clarity=1;	% To be updated!!!
	otherwise
		error(sprintf('Unknown pfType=%d!', ptOpt.pfType));
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
	if ptOpt.useParabolicFit
		line(finetunedPitchIndex*[1 1], axisLimit(3:4), 'color', 'r');
		legend('PF', 'Lower bound', 'Upper bound', 'All local max', 'Within-bound local max', 'Selected peak', 'After parabolic fit');
	else
		legend('PF', 'Lower bound', 'Upper bound', 'All local max', 'Within-bound local max', 'Selected peak');
	end
	
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);