userDataGet;

currPt=get(gca, 'CurrentPoint');
x=currPt(1,1);
y=currPt(1,2);
distance=abs(x-frameTime);
[minDist, frameIndex]=min(distance);
fprintf('\nFrame index = %d\n', frameIndex);
%frame=wave((frameIndex-1)*PP.frameSize+(1:PP.frameSize));
frame=framedY(:, frameIndex);
set(bar1H, 'xdata', frameTime(frameIndex)*[1 1]);
set(bar2H, 'xdata', frameTime(frameIndex)*[1 1]);
set(lFrameH, 'xdata', ((frameIndex-1)*(PP.frameSize-PP.overlap)+1)*[1 1]/PP.fs);
set(rFrameH, 'xdata', ((frameIndex-1)*(PP.frameSize-PP.overlap)+PP.frameSize)*[1 1]/PP.fs);

figName='Frame Plot';
figH=findobj(0, 'name', figName);
if isempty(figH)
	figure('name', figName);
else
	figure(figH);
end

userDataSet;
feval(PP.frame2pitchFcn, frame, 1, PP);