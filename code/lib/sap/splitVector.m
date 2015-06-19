function [segment, threshold] = splitVector(vector, maxLength, minLength, dispOpt)
% splitVector: Split a vector (for speech file to sentence conversion)
%	Usage: [segment, threshold] = splitVector(vector, maxLength, minLength)
%			vector: Input vector
%			maxLength: maxLength of the output segments
%			minLength: minLength of the output segments
%			segment: the output segment
%			threshold: the corresponding threshold that generates the output segment

%	Roger Jang, 20030517

if nargin<1; selfdemo; return; end
if nargin<2; maxLength=length(vector)/5; end
if nargin<3; minLength=0; end
if nargin<4; dispOpt=0; end

vector=vector(:);
maxValue=max(vector);
minValue=min(vector);
if maxValue==minValue, error('Range of the vector is zero!'); end

allThreshold=linspace(minValue, maxValue, 101);
for i=1:length(allThreshold)
	vec=vector>allThreshold(i);
	vec=[0; vec; 0];
	diffVec=diff(vec);
	index1=find(diffVec==1);
	index2=find(diffVec==-1);
	segmentLength=index2-index1;
	if max(segmentLength)<=maxLength
		threshold=allThreshold(i);
		break;
	end	
end

% 刪除長度太短的 segment
voidIndex=find(segmentLength<minLength);
index1(voidIndex)=[];
index2(voidIndex)=[];
segmentLength(voidIndex)=[];

% 傳回結果
for i=1:length(index1)
	segment(i).beginIndex=index1(i);
	segment(i).endIndex=index2(i)-1;
	segment(i).length=segmentLength(i);
end

if dispOpt
	plot(vector, '.-');
	line([1, length(vector)], threshold*[1 1], 'color', 'r');
	for i=1:length(segment)
		line([segment(i).beginIndex, segment(i).endIndex], threshold*[1 1], 'color', 'k', 'linewidth', 3); 
	end
end

% ====== Self demo
function selfdemo
t=linspace(0,7);
y=sin(t).*exp(-t/50);
y=y.^2;
maxLength=30;
minLength=5;
[segment, threshold]=feval(mfilename, y, maxLength, minLength, 1);
fprintf('maxLength = %g\n', maxLength);
fprintf('minLength = %g\n', minLength);
fprintf('%d segments identified:\n', length(segment));
for i=1:length(segment)
	fprintf('Length of segment(%d) is %g\n', i, segment(i).length);
end