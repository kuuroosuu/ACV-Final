function [localMinIndex, localMinValue, allLocalMinIndex]=amdf2ppc(amdf, fs, nbits, ptOpt, plotOpt)
% amdf2ppc: AMDF to pitch period candidates (PPC)
%	Usage: localMinIndex=amdf2ppc(amdf, fs, ptOpt, plotOpt)
%			amdf: AMDF from integer-based wave input

%	Roger Jang, 20040524

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, nbits=8; end
if nargin<4, ptOpt=ptOptSet(fs, nbits); end
if nargin<5, plotOpt=0; end

if isempty(ptOpt), ptOpt=ptOptSet(fs, nbits); end
ptOpt.dpUseLocalOptim=1;

ptOpt.minPeriod=round(ptOpt.fs/pitch2freq(ptOpt.maxPitch));
ptOpt.maxPeriod=round(ptOpt.fs/pitch2freq(ptOpt.minPitch));

% ====== т pitch period candidate 翴
if ptOpt.dpUseLocalOptim			% Method 1: Use local minima as the PPC (primary pitch candidates)
	allLocalMinIndex=localMinMex(amdf);
	localMinIndex=allLocalMinIndex; localMinValue=amdf(localMinIndex);
	if isempty(localMinIndex), amdf2ppcPlot(amdf, localMinIndex, allLocalMinIndex, ptOpt, plotOpt); return; end
else					% Method 2: Use all possible points as the PPC!
	allLocalMinIndex=1:length(amdf);
	localMinIndex=allLocalMinIndex;
	localMinValue=amdf;
end
%fprintf('length(localMinIndex)=%d\n', length(localMinIndex));

% ====== 奔 ptOpt.minPeriod ┪ ptOpt.maxPeriod  localMinIndex
index=find(localMinIndex-1<ptOpt.minPeriod | localMinIndex-1>ptOpt.maxPeriod);
localMinIndex(index)=[];  localMinValue=amdf(localMinIndex);
if isempty(localMinIndex), amdf2ppcPlot(amdf, localMinIndex, allLocalMinIndex, ptOpt, plotOpt); return; end

% ====== 玂痙 ptOpt.ppcNum  localMinIndex
%localMinValue=amdf(localMinIndex);
%[junk, index]=sort(localMinValue);
%if length(index)>ptOpt.ppcNum
%	localMinIndex=localMinIndex(index(1:min(ptOpt.ppcNum,length(index))));
%end
%localMinValue=amdf(localMinIndex);

% 秈︽繵锣传盢琌繵AMDF
%threshold=round(max(amdf)/2);
%amdfDeduct=round(max(amdf)/16);
%amdfDeduct=0;
%for i=1:length(localMinIndex)
%	for j=1:length(localMinIndex)
%		if localMinValue(i)<threshold & localMinValue(j)<threshold & localMinIndex(j)<localMinIndex(i)
%			if abs((localMinIndex(i)-1)/4-(localMinIndex(j)-1)) <= 6/4
%				if plotOpt, fprintf('传Θ 4 繵\n'); end
%				localMinValue(j)=localMinValue(j)-amdfDeduct;
%			elseif abs((localMinIndex(i)-1)/3-(localMinIndex(j)-1)) <= 6/3
%				if plotOpt, fprintf('传Θ 3 繵\n'); end
%				localMinValue(j)=localMinValue(j)-amdfDeduct;
%			elseif abs((localMinIndex(i)-1)/2-(localMinIndex(j)-1)) <= 6/2
%				if plotOpt, fprintf('传Θ 2 繵\n'); end
%				localMinValue(j)=localMinValue(j)-amdfDeduct;
%			end
%		end
%	end
%end

amdf2ppcPlot(amdf, localMinIndex, allLocalMinIndex, ptOpt, plotOpt);

% ====== Plotting
function amdf2ppcPlot(amdf, localMinIndex, allLocalMinIndex, ptOpt, plotOpt)
if plotOpt
	plot(amdf, '.-');
	title('Pitch period candidates from AMDF');
	if ~isempty(allLocalMinIndex)
		line(allLocalMinIndex, amdf(allLocalMinIndex), 'color', 'r', 'marker', 'v', 'linestyle', 'none');
	end
	if ~isempty(localMinIndex)
		line(localMinIndex, amdf(localMinIndex), 'color', 'k', 'marker', 'square', 'linestyle', 'none');
	end
	ylim=get(gca, 'ylim');
	line(ptOpt.minPeriod*[1 1], ylim, 'color', 'r');
	line(ptOpt.maxPeriod*[1 1], ylim, 'color', 'r');
	set(gca, 'xlim', [1, length(amdf)]);
end

% ====== Self demo
function selfdemo
waveFile='mountain.wav';
waveFile='Since his wife died, his life has had no meaning_70.wav';
[y, fs, nbits]=wavReadInt(waveFile);
ptOpt=ptOptSet(fs, nbits);
framedY=buffer(y, ptOpt.frameSize, ptOpt.overlap);
frameIndex=100;
frame=framedY(:, frameIndex);
frame2=frameFlip(frame);
%frame2=localAverage(frame2);
%amdf=round(frame2amdf(frame2, 128, 2, 1));	% Mimic MCU
amdf=frame2amdf(frame2, ptOpt.frameSize/2, 3);
plotOpt=1;
subplot(2,1,1);
feval(mfilename, amdf, fs, nbits, ptOpt, plotOpt);

subplot(2,1,2);
n=10;
fs=n*fs;
ptOpt=ptOptSet(fs, nbits);
newAmdf=interp1(1:length(amdf), amdf, linspace(1, length(amdf), n*length(amdf)), 'spline');
feval(mfilename, newAmdf, fs, nbits, ptOpt, plotOpt);
title('After spline interpolation of AMDF vector');
xlabel(sprintf('Frame %d of "%s"', frameIndex, waveFile));