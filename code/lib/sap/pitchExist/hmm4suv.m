function [optPath, stateProb, hmmTable]=hmm4suv(feaMat, gmmSet, transProb, plotOpt)
% hmm4suv: HMM for SU/V detection
%	Usage: [optPath, stateProb]=hmm4suv(feaMat, gmmSet, transProb, plotOpt)

%	Roger Jang, 20090928

if nargin<1, selfdemo; return; end
if nargin<4, plotOpt=1; end

[feaDim, frameNum]=size(feaMat);
gmmNum=length(gmmSet.gmm);
stateProb=zeros(frameNum, gmmNum);

for i=1:gmmNum
	[stateProb(:,i), eachProb1]=gmmEval(feaMat, gmmSet.gmm(i).gmmParam);
end

priors=[238143, 662357];	% This should be moved out as an input parameter!
logPriorProb=log(priors/sum(priors));
for i=1:frameNum
	stateProb(i,:)=stateProb(i,:)+logPriorProb(i);
end

hmmTable=zeros(frameNum, 2);
prevPos=zeros(frameNum, 2);
optPath=zeros(frameNum, 1);

% ====== Fill hmmTable
hmmTable(1,1)=stateProb(1,1);		% Initial prob for SU
hmmTable(1,2)=stateProb(1,2)-inf;	% Initial prob for V
for i=2:frameNum
	for j=1:2
		prob(1)=hmmTable(i-1,1)+transProb(1,j);
		prob(2)=hmmTable(i-1,2)+transProb(2,j);
		[maxProb, index]=max(prob);
		prevPos(i,j)=index;
		hmmTable(i,j)=stateProb(i,j)+maxProb;
	end
end

% ===== Backtrack to find the optimal path
[maxProb, optPath(end)]=max(hmmTable(end, :));
for j=frameNum-1:-1:1
	optPath(j)=prevPos(j+1, optPath(j+1));
end

if plotOpt
	subplot(5,1,1); imagesc(feaMat); axis xy
	subplot(5,1,2); plot((1:frameNum)', stateProb); legend('SU', 'V'); set(gca, 'xlim', [-inf inf]);
	subplot(5,1,3); plot(1:frameNum, optPath, 'k.-'); legend('Predcited'); set(gca, 'xlim', [-inf inf]);
	subplot(5,1,4); imagesc(stateProb'); axis xy
	subplot(5,1,5); imagesc(hmmTable'); axis xy
end
