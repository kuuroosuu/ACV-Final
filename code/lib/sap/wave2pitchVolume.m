function [pitch, volume, PP] = wave2pitchVolume(wave, duration, plotOpt, PP);
% wave2pitchVolume Get pitch and volume from a wave signal or  file

if nargin<1, selfdemo; return; end
if nargin<2, duration=5; end
if nargin<3, plotOpt=0; end
if nargin<4,
	[PP, CP]=setParam;
	PP.frame2pitchFcn='frame2pitch';
end

if ischar(wave)
	[y, PP.fs, nbits]=waveFileRead(wave);
else
	y=wave;
	PP.fs=8000;
end
y=y(1:min(duration*PP.fs, length(y)));			% 最多只取 duration 秒
framedY=buffer2(y, PP.frameSize, PP.overlap);
frameNum=size(framedY, 2);

pitch=zeros(frameNum, 1);
volume=zeros(frameNum, 1);
for i=1:frameNum,
	pitch(i)=feval(PP.frame2pitchFcn, framedY(:,i));
	volume(i)=frame2volume(framedY(:,i));
end

if plotOpt
	plotNum=3;
	subplot(plotNum,1,1);
	plot((1:length(y))/PP.fs, y);
	if ischar(wave), title(['Wave file = "', wave, '"']); end
	axis tight; grid on
	
	frameTime=(1:frameNum)*(PP.frameSize-PP.overlap)/PP.fs;
	subplot(plotNum,1,2);
	plot(frameTime, volume, '.-r');
	title('Sum of abs. magnitude');
	axis tight; grid on
	
	subplot(plotNum,1,3);
	temp=pitch;
	temp(temp==0)=nan;
	plot(frameTime, temp, '.-');
	title('Raw pitch');
	axis tight; grid on
end

% ====== Self demo
function selfdemo
waveFile='10 little indians.wav';
duration=8;
plotOpt=1;
[pitch, volume, PP]=feval(mfilename, waveFile, duration, plotOpt);
fprintf('Hit return to play the pitch:\n');
pause;
playmid(pitch/10, (PP.frameSize-PP.overlap)/PP.fs);