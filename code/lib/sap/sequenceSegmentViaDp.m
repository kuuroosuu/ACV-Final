function [normalizedMinDist, dtwPath, dtwTable, note]=sequenceSegmentViaDp(vec1, opt, showPlot)
%SequenceSegmentViaDp: Sequence segmentation via DP
%
%	Usage:
%		[minDist, dtwPath, dtwTable, note]=sequenceSegmentViaDp(vec, opt, showPlot)
%
%	Example:
%		waveFile='twinkle_twinkle_little_star.wav';
%		wObj=waveFile2obj(waveFile);
%		pfType=1;	% 0 for AMDF, 1 for ACF
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
%		ptOpt.mainFun='maxPickingOverPf';
%		pv=pitchTracking(wObj, ptOpt);
%		segment=segmentFind(pv);
%		vec=pv(segment(2).begin:segment(2).end);
%		showPlot=1;
%		opt.frameRate=wObj.fs/(ptOpt.frameSize-ptOpt.overlap);
%		opt.segmentNum=3;
%		opt.minNoteDuration=0.1;
%		[minDist, dtwPath, dtwTable, note]=sequenceSegmentViaDp(vec, opt, showPlot);
%		figure; plot((1:length(vec))/opt.frameRate, vec, '.-g');	% Plot the pv and the segmented notes
%		hold on; notePlot(note, 1, 'b', 1/opt.frameRate); hold off
%
%	Category: Note segmentation
%	Roger Jang, 20121101

if nargin<1, selfdemo; return; end
if nargin<2, opt.segmentNum=3; end
if nargin<3, showPlot=0; end

vec2=1:opt.segmentNum;
size1=size(vec1, 2);
size2=size(vec2, 2);

% ====== Construct DTW table
dtwTable=inf*ones(size1,size2);
% ====== Construct prevPos table for back tracking the optimum path

for i=1:size1
	for j=1:size2
		if nargout>1 || showPlot
			prevPos(i,j).i=-1;
			prevPos(i,j).j=-1;
		end
		prevPos(i,j).startI=-1;
	end
end

% ====== Construct the first element of the DTW table
dtwTable(1,1)=0; prevPos(1,1).startI=1;
% ====== Construct the first column of the DTW table
for j=2:size2
	dtwTable(1,j)=inf; prevPos(1,j).startI=-1;
end
% ====== Construct the first row of the DTW table (xy view)
for i=2:size1
	theVec=vec1(1:i);
	dtwTable(i,1)=sum(abs(theVec-median(theVec)));
	if nargout>1 || showPlot
		prevPos(i,1).i=i-1;
		prevPos(i,1).j=1;
	end
	prevPos(i,1).startI=1;
end
% ====== Construct all the other rows of DTW table
for i=2:size1
	for j=2:size2
		startI=prevPos(i-1,j).startI;
		if startI==-1
			dist1=inf;
		elseif startI==1
			theVec=vec1(1:i);
			dist1=sum(abs(theVec-median(theVec)));	% 0-degree path
		else
			theVec=vec1(startI:i);
			dist1=dtwTable(startI-1,j)+sum(abs(theVec-median(theVec)));	% 0-degree path
		end
		dist2=dtwTable(i-1, j-1);	% 45-degree path
		if ((i-prevPos(i-1,j-1).startI)/opt.frameRate<opt.minNoteDuration) | (dist2>=dist1)	% ====== Take 0-degree path
			dtwTable(i,j)=dist1;
			if nargout>1 || showPlot
				prevPos(i,j).i=i-1;
				prevPos(i,j).j=j;
				prevPos(i,j).startI=prevPos(i-1,j).startI;
			end
		else			% ====== Take 45-degree path
			dtwTable(i,j)=dist2;
			if nargout>1 || showPlot
				prevPos(i,j).i=i-1;
				prevPos(i,j).j=j-1;
				prevPos(i,j).startI=i;
			end
		end
	end
end

% ====== Find the last point and the min. dist
besti=size1;
minDist=dtwTable(end, end);
normalizedMinDist=minDist/length(vec1);
bestj=size2;

if nargout>1 || showPlot	% Return the optimum path
	% ====== Back track to find all the other points
	dtwPath=[besti; bestj];		% The last point in the optimum path
	nextPoint=[prevPos(dtwPath(1,1), dtwPath(2,1)).i; prevPos(dtwPath(1,1), dtwPath(2,1)).j];
	while nextPoint(1)>0 & nextPoint(2)>0
		dtwPath=[nextPoint, dtwPath];
		nextPoint=[prevPos(dtwPath(1,1), dtwPath(2,1)).i; prevPos(dtwPath(1,1), dtwPath(2,1)).j];
	end
end

for i=1:size2
	yPos=dtwPath(2,:);
	note(i).pvIndex=find(yPos==i);
	note(i).pitch=median(vec1(note(i).pvIndex));
	note(i).duration=length(note(i).pvIndex)/opt.frameRate;
	note(i).meanPitchError=mean(abs(vec1(note(i).pvIndex)-note(i).pitch));
end

% ====== Plotting if necessary
if showPlot, dtwPathPlot(vec1, vec2, dtwPath, 'auto', minDist); end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
