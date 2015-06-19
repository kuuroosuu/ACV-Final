function diffSum=pvDiffSum(pv)
% pvSimilarity: Similarity between two pitch vectors
%	Usage: [recogRate, absAveError, correctCount, pitchCount]=pvSimilarity(cPitch, tPitch, pitchTol, plotOpt)
%		cPitch: computed pitch
%		tPitch: target pitch
%		pitchTol: pitch tolerance for computing the recog. rate
%		plotOpt: 1 for plotting the result

%	Roger Jang, 20091001

segment=segmentFind(pv);
diffSum=0;
for i=1:length(segment)
	diffSum=diffSum+sum(abs(diff(pv(segment(i).begin:segment(i).end))));
end
