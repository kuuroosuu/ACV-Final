function [onset, insertCount, deleteCount]=odByVolViaFile(waveFile, odPrm, plotOpt)
% odByVolViaFile: Onset detection of tapping by volume, with file input
%	Usage: [onset, insertCount, deleteCount]=odByVolViaFile(waveFile, odPrm, plotOpt)
%
%	For example:
%		waveFile='tapping.wav';
%		odPrm=odPrmSet;
%		plotOpt=1;
%		[onset, insertCount, deleteCount]=odByVolViaFile(waveFile, odPrm, plotOpt);
%		fprintf('waveFile=%s, insertCount=%d, deleteCount=%d\n', waveFile, insertCount, deleteCount);

%	Roger Jang, 20090607

if nargin<1, selfdemo; return; end
if nargin<2, odPrm=odPrmSet; end
if nargin<3, plotOpt=1; end

[y, fs, nbits, opts, cueLabel]=wavReadInt(waveFile);
wObj=waveFile2obj(waveFile);
onset=odByVol(wObj, odPrm, plotOpt);
[insertCount, deleteCount]=insertDeleteCount(onset, cueLabel, 0.02*fs);

if plotOpt
	subplot(2,1,1);
	axisLimit=axis;
	% Display the human-transcribed cue labels
	line(cueLabel/fs, axisLimit(4)*ones(length(cueLabel),1), 'color', 'r', 'marker', 'v', 'linestyle', 'none'); 
end

% ====== Self demo
function selfdemo
waveFile='tappingNoisy.wav';
odPrm=odPrmSet;
plotOpt=1;
[onset, insertCount, deleteCount]=feval(mfilename, waveFile, odPrm, plotOpt);
fprintf('waveFile=%s, insertCount=%d, deleteCount=%d\n', waveFile, insertCount, deleteCount);
