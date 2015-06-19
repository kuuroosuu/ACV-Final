function waveData=waveDataFeaCollect(taskOpt, showPlot)
% waveDataFeaCollect: Collect all the wave files and the corresponding features
%
%	Usage:
%		waveData=waveDataFeaCollect(taskOpt, showPlot)
%
%	Example:
%		taskOpt=vdOptSet;
%		waveData=waveDataFeaCollect(taskOpt);

if nargin<1, selfdemo; return; end
if nargin<2, showPlot=0; end

waveData=recursiveFileList(taskOpt.waveDir, 'wav');
waveNum=length(waveData);
%% ====== Get all trainning data
myTic=tic;
showPlot=0;
for i=1:waveNum
	wObj=waveFile2obj(waveData(i).path);
	fprintf('%d/%d: file=%s, duration=%g sec, time=', i, waveNum, waveData(i).path, length(wObj.signal)/wObj.fs);
	theTic=tic;
	% === Collect data for each wave file
	[waveData(i).feature, waveData(i).tOutput, waveData(i).other]=taskOpt.feaExtractFun(wObj, taskOpt, showPlot);
	fprintf('%g sec\n', toc(theTic));
	if showPlot
		fprintf('Press any key to continue...'); pause; close all; fprintf('\n');
	end
end

deleteIndex=[];
for i=1:waveNum
	if isempty(waveData(i).tOutput)
		deleteIndex=[deleteIndex, i];
	end
end
fprintf('Deleting %d wave file (due to no tOutput)...\n', length(deleteIndex));
waveData(deleteIndex)=[];
fprintf('Total time=%g sec\n', toc(myTic));

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);