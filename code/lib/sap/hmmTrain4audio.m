function [hmmModel, overallRr, eachRr]=hmmTrain4audio(waveData, taskOpt, showPlot)
%hmmTrain4audio: Train HMM for audio applications
%
%	Example:
%		load vdWaveData.mat
%		taskOpt=vdOptSet;
%		hmmModel=hmmTrain4audio(waveData, taskOpt, 1);
%
%	See also hmmEval4audio.

%	Category: HMM training for audio applications
%	Roger Jang, 20130107

if nargin<1, selfdemo; return; end
if nargin<3, showPlot=0; end

%% === Transition probability
%tic, fprintf('Computing transition probability...');
tOutput=[waveData.tOutput];
feature=[waveData.feature];
tempOutput=tOutput;	% For computing transition prob.
allTrans = enframe(tOutput, 2, 1);
if isfield(waveData(1).other, 'invalidIndex')
	temp=[waveData.other]; invalidIndex=[temp.invalidIndex];
	feature(:, invalidIndex)=[];
	tOutput(:, invalidIndex)=[];
	tempOutput(invalidIndex)=0;		% Break the sequence for computing trans. prob.
	validIndex=all(allTrans);	% Keep valid columns which contain no zeros
	allTrans=allTrans(:, validIndex);
end
classNum=length(unique(tOutput));
transNum=size(allTrans,2);
[uniqueBigram, bigramCount]=rowCount(allTrans');
countMat=zeros(classNum);
for i=1:size(uniqueBigram,1)
	countMat(uniqueBigram(i,1), uniqueBigram(i,2))=bigramCount(i);
end
hmmModel.transProb=bsxfun(@times, countMat, 1./sum(countMat,2));	% Transition probability
hmmModel.transProb=log(hmmModel.transProb);	% Log-based transition prob. (How come this line reduces the accuracy???)
%fprintf('time=%g sec\n', toc);
%% === State probability
%tic, fprintf('Building GMM models...');

opt=gmmTrain('defaultOpt');
opt.config.gaussianNum=taskOpt.gaussianNum;
for i=1:classNum
	theFeature=feature(:, tOutput==i);
	[hmmModel.gmm{i}, logLike{i}] = gmmTrain(theFeature, opt);
end
%fprintf('time=%g sec\n', toc);

% === Compute RR for each file
waveNum=length(waveData);
waveData(1).cOutput=[]; waveData(1).rr=[];	% Add the fields to be attached
for i=1:waveNum
	waveData(i)=hmmEval4audio(waveData(i), taskOpt, hmmModel);
end
eachRr=[waveData.rr];
cOutput=[waveData.cOutput];
tOutput=[waveData.tOutput];
overallRr=sum(cOutput==tOutput)/length(tOutput);
	
if showPlot
	subplot(211);
	bar(eachRr*100); grid on
	ylabel('Recog. rate (%)');
	title(sprintf('Overall inside-test accuracy=%g%% for files in %s', overallRr*100, strPurify4label(taskOpt.waveDir)));
	if waveNum<=20
		set(gca, 'xtick', 1:waveNum);
		xTickLabelRename({waveData.name});
		xTickLabelRotate(330, 10, 'left');
	else
		xlabel('Wave file indices');
	end
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
