function pitchLabel(au, ptOpt)
% pitchLabel: Label the pitch of a given audio file
%
%	Usage:
%		pitchLabel(waveFile, ptOpt)
%
%	Example:
%		waveFile='twinkleTwinkleLittleStar.wav';
%		pitchLabel(waveFile);

%	Roger Jang 20040911, 20100920, 20150409

if nargin<1, selfdemo; return; end
if isstr(au), au=myAudioRead(au); end
if nargin<2
	ptOpt=ptOptSet(au.fs, au.nbits);
	ptOpt.pfLen=ptOpt.frameSize;
end
ptOpt.waveFile=au.file;
ptOpt.targetPitchFile=[ptOpt.waveFile(1:end-3), 'pv'];
ptOpt.frame2pitchFcn='frame2pitch4labeling';

% ====== Apply low-pass filter
%cutOffFreq=500;	% Cutoff frequency
%filterOrder=3;		% Order of filter
%[b, a]=butter(filterOrder, cutOffFreq/(fs/2), 'low');
%y=filter(b, a, y);
%%y=meanFilter(y, 9);

wave2pitch(au, ptOpt, 1);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
strEval(mObj.example);