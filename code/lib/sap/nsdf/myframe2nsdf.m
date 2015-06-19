function out=myframe2sndf(frame)

maxShift = length(frame);
frameSize=length(frame);
out=zeros(maxShift, 1);

% moving base = whole frame, but normalized by the overlap area
for i=1:maxShift
	tempAcf = dot(frame(1:frameSize-i+1), frame(i:frameSize));
	tempSDF = frame(1:frameSize-i+1)'*frame(1:frameSize-i+1) + frame(i:frameSize)'*frame(i:frameSize)+0.1;
	out(i) = 2*tempAcf/tempSDF;		% normalization
end