function [pitch, clarity, pitchIndex, pf]=frame2pitch(frame, ptOpt, plotOpt)
% frame2pitch: Frame to pitch conversion using PF
%	Usage: [pitch, clarity, pitchIndex]=frame2pitch(frame, ptOpt, plotOpt)
%
%	See also frame2pf, frame2acf, frame2amdf.

%	Roger Jang 20070209, 20100617

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

frame=frameZeroMean(frame, 6);
frameSize=length(frame);
switch (ptOpt.pfType)
	case 0	% AMDF
		method=2;
		pf=frame2amdfMex(frame, frameSize, method);
	%	pf=pf.*(0:length(pf)-1)';
		[pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt);
	case 1	% ACF
		method=2;
		pf=frame2acfMex(frame, frameSize, method);
	%	if method==1
	%		pfWeight=1+linspace(0, ptOpt.alpha, length(pf))';
	%		pf=pf.*pfWeight;	% To avoid double pitch error (esp for violin). 20110416
	%	end
		if method==2
			pfWeight=1-linspace(0, ptOpt.alpha, length(pf))';	% alpha is less than 1.
			pf=pf.*pfWeight;	% To avoid double pitch error (esp for violin). 20110416
		end
		[pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt);
	case 2	% NSDF
		method=1;
		pf=frame2nsdfMex(frame, frameSize, method);
		[pitch, clarity, pitchIndex]=pf2pitch(pf, ptOpt);
	case 3	% HPS
		zeroPaddedFrameSize=16*length(frameSize);
		[pf, freq]=frame2hps(frame, ptOpt.fs, zeroPaddedFrameSize);
		[maxValue, maxIndex]=max(pf);
		pitch=freq2pitch(freq(maxIndex));
		clarity=0.8;
		pitchIndex=nan;
	otherwise
		error('Unknown pfType!');
end

if plotOpt
	subplot(2,1,1); plot(frame); axis tight; title('Frame');
	subplot(2,1,2); [pitch, clarity, ppcIndex]=pf2pitch(pf, ptOpt, plotOpt); title(sprintf('PF (pfType=%d)', ptOpt.pfType));
end

% ====== Self demo
function selfdemo
waveFile='greenOil.wav';
[y, fs, nbits]=wavReadInt(waveFile);
frameMat=buffer2(y, 256, 0);
frame=frameMat(:, 250);
pfType=2;
ptOpt=ptOptSet(fs, nbits, pfType);
plotOpt=1;
[pitch, clarity, pitchIndex]=frame2pitch(frame, ptOpt, plotOpt);
