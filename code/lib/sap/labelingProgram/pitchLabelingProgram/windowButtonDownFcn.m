userDataGet;

currPt=get(gca, 'CurrentPoint');
x=currPt(1,1);
y=currPt(1,2);
distance=abs(x-frameTime);
[minDist, frameIndex]=min(distance);
%fprintf('\nFrame index = %d\n', frameIndex);
%frame=wave((frameIndex-1)*ptOpt.frameSize+(1:ptOpt.frameSize));
frame=frameMat(:, frameIndex);
set(bar1H, 'xdata', frameTime(frameIndex)*[1 1]);
set(bar2H, 'xdata', frameTime(frameIndex)*[1 1]);
set(lFrameH, 'xdata', ((frameIndex-1)*(ptOpt.frameSize-ptOpt.overlap)+1)*[1 1]/ptOpt.fs);
set(rFrameH, 'xdata', ((frameIndex-1)*(ptOpt.frameSize-ptOpt.overlap)+ptOpt.frameSize)*[1 1]/ptOpt.fs);

figTagName='My Frame Plot';
figH=findobj(0, 'tag', figTagName);
if isempty(figH)
	figH=figure('tag', figTagName);
else
	figure(figH);
end
set(figH, 'name', sprintf('Current frame index=%d', frameIndex));

userDataSet;
frame2pitch4labeling(frame, ptOpt, 1);