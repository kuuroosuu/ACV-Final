function w=buffer(w, frameSize, overlap)
%BUFFER Wave to frame matrix conversion
%	Usage: w=buffer(w, frameSize, overlap)
%		w: wave object
%		frameSize: frame size
%		overlap: overlap

%	Roger Jang, 20060815

step = frameSize-overlap;
frameCount = floor((length(w.signal)-overlap)/step);
w.frameMatrix=zeros(frameSize, frameCount);
for i=1:frameCount,
	startIndex = (i-1)*step+1;
	w.frameMatrix(:, i)=w.signal(startIndex:(startIndex+frameSize-1));
end
w.frameSize=frameSize;
w.overlap=overlap;