function wObj=hmmEval4audio(wObj, taskOpt, hmmModel, showPlot)
% hmmEval4audio: HMM evaluation for vibrato detection
%
%	Usage:
%		wObj=hmmEval4audio(wObj);
%		wObj=hmmEval4audio(wObj, taskOpt, hmmModel);
%		wObj=hmmEval4audio(wObj, taskOpt, hmmModel, showPlot);
%
%	Description:
%		The output wObj will be attached with three fields: tOutput, cOutput, and rr.
%
%	Example:
%		waveFile='D:\dataset\vibrato\female\combined-female.wav';
%		taskOpt=vdOptSet;
%		load vdHmmModel.mat	% Obtain hmmModel
%		wObj=hmmEval4audio(waveFile, taskOpt, hmmModel, 1);

%	Category: HMM evaluation for audio applications
%	Roger Jang, 20130106

if nargin<1, selfdemo; return; end
if nargin<4, showPlot=0; end

if ischar(wObj), wObj=waveFile2obj(wObj); end
%% Feature extraction
if ~isfield(wObj, 'feature')
	[wObj.feature, wObj.tOutput, wObj.other]=taskOpt.feaExtractFun(wObj, taskOpt);
end
%% HMM evaluation
classNum=length(hmmModel.gmm);
dataNum=size(wObj.feature,2);
map=zeros(classNum, dataNum);
for i=1:classNum
	map(i, :)=gmmEval(wObj.feature, hmmModel.gmm{i});
end
[optimValue, dpPath, dpTable]=dpOverMap(map, hmmModel.transProb);
wObj.cOutput=dpPath(2,:);
%if ~isempty(wObj.other.cueTime)
%	interval=rangeSearchBin(wObj.other.sFrameTime, wObj.other.cueTime);
%	wObj.tOutput=2-mod(interval, 2);	% 1: non-vibrato, 2: vibrato
%end
wObj.rr=sum(wObj.tOutput==wObj.cOutput)/length(wObj.tOutput);
%% Plotting
if showPlot
	taskOpt.hmmPlotFun(wObj, taskOpt);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
