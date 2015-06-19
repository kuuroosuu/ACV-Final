function [outsideRr, cvData]=hmmPerfLoo4audio(waveData, taskOpt, showPlot)
%perfLoo: Performance evaluation via LOO
%
%	Usage:
%		[outsideRr, cvData]=perfLoo(waveData)
%
%	Example:
%		load vdWaveData.mat
%		[outsideRr, cvData]=hmmPerfLoo4audio(waveData, 1);

%	Category: Performance evaluation
%	Roger Jang, 20130107

if nargin<1, selfdemo; return; end
if nargin<3, showPlot=0; end

waveNum=length(waveData);
for i=1:waveNum
	if showPlot, fprintf('%d/%d: file=%s\n', i, waveNum, waveData(i).path); end
	loopTic=tic;
	theWaveData=waveData;
	theWaveData(i)=[];
	[cvData(i).hmmModel, cvData(i).insideRr, eachRr]=hmmTrain4audio(theWaveData, taskOpt);
	temp=hmmEval4audio(waveData(i), taskOpt, cvData(i).hmmModel);
	cvData(i).cOutput=temp.cOutput; cvData(i).tOutput=temp.tOutput;
	cvData(i).outsideRr=sum(cvData(i).cOutput==cvData(i).tOutput)/length(cvData(i).tOutput);
	if showPlot, fprintf('\toutsideRr=%g%%, time=%g sec\n', cvData(i).outsideRr*100, toc(loopTic)); end
end
cOutput=[cvData.cOutput];
tOutput=[cvData.tOutput];
outsideRr=sum(cOutput==tOutput)/length(tOutput);
if showPlot, fprintf('Overall LOO accuracy=%g%%\n', outsideRr*100); end

if showPlot
%	subplot(211);
	thisPlot(cvData, waveData, outsideRr, taskOpt);
%	[sorted, index]=sort([cvData.outsideRr]); cvData2=cvData(index); waveData2=waveData(index);
%	subplot(212);
%	thisPlot(cvData2, waveData2, outsideRr, taskOpt);
end

% ====== Plot for this function
function thisPlot(cvData, waveData, outsideRr, taskOpt)
plot(1:length(cvData), [cvData.insideRr]*100, '.-', 1:length(cvData), [cvData.outsideRr]*100, '.-'); grid on
legend('Inside-test RR', 'Outside-test RR');
ylabel('Recog. rate (%)');
title(sprintf('Overall LOO accuracy=%g%% for files in %s', outsideRr*100, strPurify4label(taskOpt.waveDir)));
set(gca, 'xlim', [1, length(cvData)]);
set(gca, 'xtick', 1:length(cvData));
if length(cvData)<=20
	xTickLabelRename({waveData.name});
	xTickLabelRotate(330, 10, 'left');
else
	xlabel('File indices');
	labelOfIndex=cellfun(@num2str, num2cell(1:length(cvData)), 'uniformoutput', false);
	xTickLabelRename(labelOfIndex);
end
	
% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
