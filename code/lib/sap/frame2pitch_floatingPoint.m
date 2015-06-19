function [pitch, frame2, frame3, amdf]=frame2pitch(frame, plotOpt, PP, mainWindow, correctPitch);
%frame2pitch2: Frame to pitch conversion (Same as frame2pitch, except that everything is floating-point)
%	Usage: pitch=frame2pitch(frame, plotOpt, PP); 
%		frame: Each element is unsigned integer between 0 and 255 (inclusive).
%		plotOpt: 1 for plotting, 0 for not plotting
%		pitch: Output pitch in semitone

%	Roger Jang, 20021201

if nargin<1, selfdemo; return; end
if nargin<2, plotOpt=0; end
if nargin<3, [PP, CP]=setParam; end
if nargin<4, mainWindow=0; end
if nargin<5, correctPitch=0; end

%saveAscii(frame, 'text/frame.txt')

PP.maxFreq=semitone2freq(PP.maxPitch);
PP.minFreq=semitone2freq(PP.minPitch);
fs=PP.fs;
maxShift=PP.maxShift;
frameSize=length(frame);

% 如果音量太小，則直接回傳 0
%if sum(abs(frame-128)) < 9*256,
%	pitch = 0;
%	return;
%end

% 求 frame 的平均值為
average=mean(frame);
frame=frame-average;

% frame2: 如果前半音框音量小於後半音框，進行鏡射
frame2=frameFlip(frame, plotOpt);

% frame3: 平滑化
frame3=localAverage(frame2);	% Low-pass filter implemented as local average

% ====== 計算 AMDF 曲線
amdf=frame2amdf(frame3, maxShift);

% ====== Find ROI
beginIndex=ceil(PP.fs/PP.maxFreq);
endIndex=min(floor(PP.fs/PP.minFreq), PP.maxShift);
% Move beginIndex forward until amdf is going down
for i=beginIndex:endIndex-1,
	if amdf(i)>amdf(i+1),
		break;
	end
end
beginIndex=i+1;
if beginIndex==endIndex,
%	fprintf('beginIndex==endIndex ==> pitch = 0！\n');
	pitch=0;
	plotFrame(plotOpt, frame, frame3, amdf, [], [], beginIndex, endIndex, mainWindow, correctPitch, PP.maxShift);
	return;
end

% ====== Find local minima in ROI
localMinIndex=[];
for i=beginIndex+1:endIndex-1
	if amdf(i-1)>amdf(i) & amdf(i)<=amdf(i+1)
		localMinIndex=[localMinIndex, i];	
	end
end

n=20;
n=100;
if length(localMinIndex)>=n,	% 應該是氣音（但對女聲會誤判）
	pitch=0;
	if plotOpt, fprintf('Local minima 超過 %d 個 ==> pitch = 0！\n', n); end
	plotFrame(plotOpt, frame, frame3, amdf, localMinIndex, [], beginIndex, endIndex, mainWindow, correctPitch, PP.maxShift);
	return
end

% ====== 找最小值 in ROI
roi=amdf(beginIndex:endIndex);	% region of interest
[minValue, minIndex]=min(roi);
minIndex=minIndex+beginIndex-1;


% ===== 從 minIndex 往回找可能出現的 2, 3, 4, 5, 6 倍頻
[maxValue, maxIndex]=max(roi);
difthreshold=minValue+floor((maxValue-minValue)/8);
for i=1:length(localMinIndex)
	if amdf(localMinIndex(i)) <= difthreshold
		if abs((minIndex-1)/2-(localMinIndex(i)-1)) <= 6/2
			if plotOpt, fprintf('代換成 2 倍頻！\n'); end
			minIndex=localMinIndex(i);
			break;
		elseif abs((minIndex-1)/3-(localMinIndex(i)-1)) <= 6/3
			if plotOpt, fprintf('代換成 3 倍頻！\n'); end
			minIndex=localMinIndex(i);
			break;
		elseif abs((minIndex-1)/4-(localMinIndex(i)-1)) <= 6/4
			if plotOpt, fprintf('代換成 4 倍頻！\n'); end
			minIndex=localMinIndex(i);
			break;
		elseif abs((minIndex-1)/5-(localMinIndex(i)-1)) <= 6/5
			if plotOpt, fprintf('代換成 5 倍頻！\n'); end
			minIndex=localMinIndex(i);
			break;
		elseif abs((minIndex-1)/6-(localMinIndex(i)-1)) <= 6/6
			if plotOpt, fprintf('代換成 6 倍頻！\n'); end
			minIndex=localMinIndex(i);
			break;
		end
	end
end


freq=PP.fs/(minIndex-1);	% 取第一個 minimum 來計算基頻
pitch = 69+12*log2(freq/440);
% ====== Plot all related information
plotFrame(plotOpt, frame, frame3, amdf, localMinIndex, minIndex, beginIndex, endIndex, mainWindow, correctPitch, PP.maxShift);

% ====== Plot related information
function plotFrame(plotOpt, frame, frame3, amdf, localMinIndex, minIndex, beginIndex, endIndex, mainWindow, correctPitch, maxShift)
if plotOpt,
	frameSize=length(frame);
	subplot(2,1,1);
%	plot(1:length(frame), frame, '.-', 1:length(frame3), frame3, '.-'); grid on
	frameH=plot(1:length(frame), frame, '.-'); grid on
	set(frameH, 'tag', 'frame');
	title(['Frame (frameSize=', int2str(frameSize), ')']);
	set(gca, 'xlim', [1 frameSize]);
	subplot(2,1,2);
	amdfH=plot(amdf, '.-'); grid on
	set(amdfH, 'tag', 'amdf');
	amdfAxisH=gca;
	set(gca, 'tag', 'amdfAxis');
	title(['amdf (maxShift=', int2str(maxShift), ')']);
	set(gca, 'xlim', [1 maxShift]);
	if ~isempty(localMinIndex)
		localMinIndexH=line(localMinIndex, amdf(localMinIndex), 'linestyle', 'none', 'color', 'k', 'marker', 'o', 'tag', 'localMinIndexH');
	end
	if ~isempty(minIndex)
		minIndexH=line(minIndex, amdf(minIndex), 'linestyle', 'none', 'color', 'r', 'marker', 'o', 'tag', 'minIndexH');
		if correctPitch==0
			minIndex=nan;
		else
			minIndex=round(1+8000/(440*2^((correctPitch/10-69)/12)));
		end
		manualBarH=line(minIndex*[1 1], get(amdfAxisH, 'ylim'), 'color', 'm', 'erase', 'xor', 'tag', 'manualBar');
	end
	line(beginIndex*[1 1], get(gca, 'ylim'), 'color', 'r');
	line(  endIndex*[1 1], get(gca, 'ylim'), 'color', 'r');
	set(gcf, 'name', mfilename);
	
	% 設定滑鼠按鈕的反應動作，以便讓使用者修正 pitch
	if mainWindow>0;
		set(gcf, 'userdata', mainWindow);
		frameWinMouseAction;
	end
end


% ====== selfdemo
function selfdemo
waveFile='兩隻老虎.wav';
[y, PP.fs, nbits]=wavread(waveFile);
y=y*(2^nbits/2);

frameSize=256;
overlap=-256;
framedY=buffer2(y, frameSize, overlap);
frame=framedY(:, 100);
plotOpt=1;
feval(mfilename, frame,  plotOpt);