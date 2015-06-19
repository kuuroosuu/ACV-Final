function [pitchRange, pitchDiffMax]=pitchStats(pitch)

segment=segmentFind(pitch);
for j=1:length(segment)
	segment(j).maxPitchDiff=abs(diff(pitch(segment(j).begin:segment(j).end)));
end
pitchDiffMax=max([segment.maxPitchDiff]);
temp=pitch;
temp(temp==0)=nan;
pitchRange=max(temp)-min(temp);