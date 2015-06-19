function sampleRange=frame2sampleRange(frameRange, frameSize, overlap)

start = (frameRange(1)-1)*(frameSize-overlap)+1;
stop  = (frameRange(2)-1)*(frameSize-overlap)+frameSize;
sampleRange = [start, stop];