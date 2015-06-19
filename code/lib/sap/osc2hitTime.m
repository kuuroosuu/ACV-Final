function cHitTime=osc2hitTime(osc, opt, showPlot)
% osc2hitTime: Compute hit time from OSC (onset strength curve)
%
%	Usage:
%		cHitTime=osc2hitTime(osc, opt)
%
%	Example:
%		waveFile='D:\dataSet\igs-鈊象電子\1.歌曲\變拍\560_愛你愛到死(也有變bpm)_74.000.wav';
%		waveFile='song01s5.wav';
%		wObj=waveFile2obj(waveFile);
%		wObj.signal=mean(wObj.signal, 2);
%		oscOpt=wave2osc('defaultOpt');
%		osc=wave2osc(wObj, oscOpt);
%		hitTimeOpt=osc2hitTime('defaultOpt');
%		cHitTime=osc2hitTime(osc, hitTimeOpt, 1);

%	Roger Jang, 20120605, 20130707

if nargin<1, selfdemo; return; end
if nargin<3, showPlot=0; end

if nargin<1, selfdemo; return; end
% ====== Set the default options
if nargin==1 && ischar(osc) && strcmpi(osc, 'defaultOpt')
	cHitTime.minOscRatio=0.04;		% Was 0.069
	cHitTime.minHitSeparation=0.139;
	return
end
if nargin<2|isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

novelty=osc.signal;
noveltyTh=max(novelty)*opt.minOscRatio;
onset1=novelty>noveltyTh & localMax(novelty);

onset2=onset1;
halfWinWidth=round(opt.minHitSeparation/diff(osc.time(1:2)));
index=find(onset2);
for i=index
	startIndex=max(i-round(halfWinWidth), 1);
	endIndex=min(i+round(halfWinWidth), length(novelty));
	[junk, maxIndex]=max(novelty(startIndex:endIndex));
	if maxIndex+startIndex-1~=i
		onset2(i)=0;
	end
end
cHitTime=osc.time(onset2);

if showPlot
	plot(osc.time, novelty, '.-'); set(gca, 'xlim', [-inf inf]);
	xlabel('Time (sec)'); ylabel('Volume');
	line([osc.time(1), osc.time(end)], noveltyTh*[1 1], 'color', 'r');
	line(osc.time(onset1), novelty(onset1), 'marker', '.', 'color', 'g', 'linestyle', 'none');
	line(osc.time(onset2), novelty(onset2), 'marker', '^', 'color', 'k', 'linestyle', 'none');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
