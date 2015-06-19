function [recogRate, absAveError, correctCount, pitchCount]=pvSimilarity(cPitch, tPitch, pitchTol, plotOpt)
% pvSimilarity: Similarity between two pitch vectors
%	Usage: [recogRate, absAveError, correctCount, pitchCount]=pvSimilarity(cPitch, tPitch, pitchTol, plotOpt)
%		cPitch: computed pitch
%		tPitch: target pitch
%		pitchTol: pitch tolerance for computing the recog. rate
%		plotOpt: 1 for plotting the result

%	Roger Jang, 20091001

if nargin<1, selfdemo; return; end
if nargin<4, plotOpt=0; end

cPitch=cPitch(:);
tPitch=tPitch(:);
% Take care of length difference
if length(cPitch)~=length(tPitch)
	fprintf('Warning: Pitch lengths mismatch in %s!\n', mfilename);
end
len=min(length(cPitch), length(tPitch));
cPitch=cPitch(1:len);
tPitch=tPitch(1:len);

absDiff=abs(cPitch-tPitch);
nonZeroIndex=(cPitch~=0) & (tPitch~=0);
outOfBoundIndex=absDiff>pitchTol;
errorIndex=nonZeroIndex & outOfBoundIndex;
errorCount=sum(errorIndex);
absAveError=mean(absDiff(nonZeroIndex));
pitchCount=sum(nonZeroIndex);
correctCount=pitchCount-errorCount;
recogRate=correctCount/pitchCount;

if plotOpt
	pitchNum=length(tPitch);
	tPitch(tPitch==0)=nan;
	plot(1:pitchNum, tPitch, 'o-g', 1:pitchNum, cPitch, '.-b');
	index=find(errorIndex);
	line(index, cPitch(index), 'marker', 'x', 'linestyle', 'none', 'color', 'r');
	grid on;
	legend('Target pitch', 'Computed pitch', 'Erroneous pitch', 'location', 'northOutside', 'orientation', 'horizontal');
	xlabel(sprintf('Raw pitch accuracy=%g%%, pitchTol=%g, correctCount=%d, pitchCount=%d', recogRate*100, pitchTol, correctCount, pitchCount));
end

% ====== Self demo
function selfdemo
cPitch=[97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 56.9352 56.6997 56.9352 57.174 58.6804 58.1635 97.5248 97.5248 97.5248 97.5248 63.4868 63.4868 62.8078 60.6214 59.7627 58.1635 58.1635 58.1635 57.6617 49.5248 49.3709 61.2184 48.6214 48.4753 48.4753 49.2184 61.5248 59.7627 72.3304 72.9173 93.174 93.174 58.42 56.238 54.3004 54.3004 54.3004 54.7126 56.0117 57.174 57.9108 60.3304 59.7627 58.6804 57.4161 56.4673 55.7883 55.5677 56.0117 56.4673 56.4673 56.6997 58.1635 58.9447 58.42 58.1635 56.6997 56.0117 56.0117 56.6997 55.5677 57.174 57.6617 57.6617 57.6617 56.9352 55.5677 55.35 55.7883 56.0117 56.0117 56.0117 56.0117 56.4673 56.4673 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248 97.5248];
tPitch=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 58.680365 56.9352 56.699654 56.9352 57.173995 58.420025 57.910801 0 0 0 0 0 63.486821 62.478049 60.917323 59.485736 58.420025 58.163541 58.163541 57.661699 0 0 61.218415 60.621378 60.330408 60.621378 61.218415 61.524836 0 0 0 0 60.917323 58.944681 55.788268 54.300404 54.300404 54.505286 54.922471 55.788268 56.9352 57.910801 59.762739 59.485736 58.420025 57.416129 56.699654 56.011656 55.567726 56.237965 56.467271 56.467271 56.9352 58.163541 58.944681 58.680365 58.163541 56.9352 56.237965 56.011656 56.467271 55.567726 0 0 0 58.163541 56.699654 55.567726 55.349958 55.788268 56.011656 55.788268 56.011656 56.011656 56.467271 56.467271 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
pitchTol=1;
plotOpt=1;
feval(mfilename, cPitch, tPitch, pitchTol, plotOpt);