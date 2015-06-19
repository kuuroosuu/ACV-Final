function pitch=wave2pitch(wave, fs, nbits, plotOpt, PP);
% wave2pitch: 由 wave 資料產生音高向量
%	Usage: pitch=wave2pitch(wave, fs, nbits, plotOpt, PP); 
%		wave: Each element is a signed integer.
%		plotOpt: 1 for plotting, 0 for not plotting
%		PP: pitch parameters
%		pitch: 電腦計算之音高

%	Roger Jang, 20021201

if nargin<1, selfdemo; return; end
if nargin<4, plotOpt=0; end
if nargin<5, PP=setPitchPrm(fs, nbits); end


wave=wave-mean(wave);
framedY=buffer2(wave, PP.frameSize, PP.overlap);
frameNum=size(framedY, 2);
pitch=zeros(frameNum, 1);
volume=zeros(frameNum, 1);

% ====== 計算音量
volume=frame2volume(framedY);
% ====== 計算音量門檻值
%PP.volTh=getVolTh(volume, PP.frameSize, 4);
PP.volTh=mean(volume)-0.8*var(volume)^0.5;  % modified by felix 20061218
% ====== 計算音高（會用到音量門檻值）
for i=1:frameNum,
	pitch(i)=feval(PP.frame2pitchFcn, framedY(:,i), 0, PP);
end
plotOpt = 0; %add by felix 20061218
% ====== 畫圖
if plotOpt
	% ====== Plot wave form
	frameTime=(1:frameNum)*(PP.frameSize-PP.overlap)/PP.fs;
	plotNum=3;
	waveAxisH=subplot(plotNum,1,1);
	waveH=plot((1:length(wave))/PP.fs, wave);
	if isfield(PP, 'waveFile')
		titleStr=['Wave file = "', PP.waveFile, '"'];
		titleStr=strrep(titleStr, '_', '\_');
		title(titleStr);
	end
	axis([1/fs, length(wave)/fs, -2^PP.nbits/2, 2^PP.nbits/2-1]); grid on
	lFrameH=line(nan*[1 1], get(waveAxisH, 'ylim'), 'erase', 'xor', 'color', 'r');
	rFrameH=line(nan*[1 1], get(waveAxisH, 'ylim'), 'erase', 'xor', 'color', 'r');
	% ====== Plot volume
	volumeAxisH=subplot(plotNum,1,2);
	plot(frameTime, volume, '.-');
	title('Sum of abs. magnitude');
	line([min(frameTime), max(frameTime)], PP.volTh*[1 1], 'color', 'r');
	set(gca, 'xlim', [-inf inf]); grid on
	bar1H=line(nan*[1 1], get(gca, 'ylim'), 'erase', 'xor', 'color', 'r');
	% ====== Plot 電腦算的 pitch & 手動標示的 pitch
	pitchAxisH=subplot(plotNum,1,3);
	temp=pitch; temp(temp==0)=nan;
	pitchH=plot(frameTime, temp, 'o-', 'color', 'g');		% 電腦算的 pitch
	titleStr='電腦計算之音高';
	targetPitch=pitch;
	if isfield(PP, 'targetPitchFile') & exist(PP.targetPitchFile)	% 正確音高所在的檔案
		titleStr='電腦計算（綠色）和手動標示（黑色）之音高';
		targetPitch=asciiRead(PP.targetPitchFile);		% 從音高檔案讀出手動標示之正確音高
	end
	if length(targetPitch)>frameNum, targetPitch=targetPitch(1:frameNum); end	% 由於 buffer (used before) 和 buffer2 (used now) 不同所造成的差異
	temp=targetPitch; temp(temp==0)=nan;
	targetPitchH=line(frameTime, temp, 'color', 'k', 'marker', '.');		% 正確音高
	title(titleStr);
	axis tight; grid on
	bar2H=line(nan*[1 1], get(gca, 'ylim'), 'erase', 'xor', 'color', 'r');
	% ====== Place the figure
	screenSize=get(0, 'screensize');
	currPos=get(gcf, 'position');
	set(gcf, 'position', [screenSize(3)-currPos(3), 30, currPos(3:4)]);
	% ====== Record user data
	userDataSet;
	set(gcf, 'WindowButtonDownFcn', 'WindowButtonDownFcn');
	uicontrol('string', 'Play Wave', 'Callback', 'userDataGet; wavplay(wave/(2^PP.nbits/2), PP.fs, ''async'')', 'position', [20, 10, 60, 20]);
	uicontrol('string', 'Play Pitch1', 'Callback', 'userDataGet; pitchPlay(pitch, PP.fs/(PP.frameSize-PP.overlap));', 'position', [100, 10, 60, 20]);
	uicontrol('string', 'Play Pitch2', 'Callback', 'userDataGet; pitchPlay(targetPitch, PP.fs/(PP.frameSize-PP.overlap));', 'position', [180, 10, 60, 20]);
%	uicontrol('string', 'Save Pitch', 'Callback', 'savePitch', 'position', [260, 10, 60, 20]);
end

% ====== self demo
function selfdemo
waveFile='倫敦鐵橋垮下來_不詳_0.wav';
[y, fs, nbits]=wavread(waveFile);
wave=y*2^nbits/2;
PP=setPitchPrm(fs, nbits);
PP.waveFile=waveFile;
PP.targetPitchFile=[waveFile(1:end-3), 'pitch'];
feval(mfilename, wave, fs, nbits, 1, PP);