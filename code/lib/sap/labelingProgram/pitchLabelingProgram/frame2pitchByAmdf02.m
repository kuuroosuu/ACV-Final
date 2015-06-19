function [pitch, frame2, amdf]=frame2pitch4labeling(frame, plotOpt, ptOpt);
% frame2pitch4labeling: 由一個音框計算一點音高
%	Usage: pitch=frame2pitch4labeling(frame, ptOpt, plotOpt); 
%		frame: Each element is unsigned integer between 0 and 255 (inclusive).
%		plotOpt: 1 for plotting, 0 for not plotting
%		pitch: Output pitch in semitone

%	Roger Jang, 20021201

if nargin<1, selfdemo; return; end
if nargin<2, plotOpt=0; end
if nargin<3, ptOpt=ptOptSet(8000); end

userDataGet;

% 求 frame 的平均值
average=mean(frame);
frame=frame-average;

% frame2: 如果前半音框音量小於後半音框，進行鏡射
frame2=frameFlip(frame, plotOpt);

% ====== 計算 AMDF 曲線
acf=frame2acf(frame2, length(frame2));
%acf=frame2acfMex(frame2);
amdf=frame2amdf(frame2, length(frame2));
%amdf=frame2amdfMex(frame2);

% ====== Find ROI (region of interest)
beginIndex=ceil(ptOpt.fs/ptOpt.maxFreq);
endIndex=min(floor(ptOpt.fs/ptOpt.minFreq), ptOpt.maxShift);

% ====== Find local minima in ROI
localMinIndex=localMinMex(amdf);
localMinCount=length(localMinIndex);
localMinIndex(localMinIndex>endIndex)=[];	% Remove out-of-bound local min. 
localMinIndex(localMinIndex<beginIndex)=[];	% Remove out-of-bound local min.
if isempty(localMinIndex)
	pitch=0;
	return;
end
localMinValue=amdf(localMinIndex);
% ===== Find the min. among all local min.
[minValue, ind]=min(localMinValue);
minIndex=localMinIndex(ind);

% ===== 從 minIndex 往回找可能出現的 2, 3, 4, 5, 6 倍頻
roi=amdf(beginIndex:min(endIndex, length(amdf)));	% region of interest
minValue=min(roi);
maxValue=max(roi);
difthreshold=minValue+(maxValue-minValue)/8;
if ptOpt.checkMultipleFreq
	for i=1:length(localMinIndex)
		if amdf(localMinIndex(i))<=difthreshold
			if abs((minIndex-1)/6-(localMinIndex(i)-1)) <= 6/6
				if plotOpt, fprintf('代換成 6 倍頻！\n'); end
				minIndex=localMinIndex(i);
				break;
			elseif abs((minIndex-1)/5-(localMinIndex(i)-1)) <= 6/5
				if plotOpt, fprintf('代換成 5 倍頻！\n'); end
				minIndex=localMinIndex(i);
				break;
			elseif abs((minIndex-1)/4-(localMinIndex(i)-1)) <= 6/4
				if plotOpt, fprintf('代換成 4 倍頻！\n'); end
				minIndex=localMinIndex(i);
				break;
			elseif abs((minIndex-1)/3-(localMinIndex(i)-1)) <= 6/3
				if plotOpt, fprintf('代換成 3 倍頻！\n'); end
				minIndex=localMinIndex(i);
				break;
			elseif abs((minIndex-1)/2-(localMinIndex(i)-1)) <= 6/2
				if plotOpt, fprintf('代換成 2 倍頻！\n'); end
				minIndex=localMinIndex(i);
				break;
			end
		end
	end
end

freq=ptOpt.fs/(minIndex-1);
pitch=freq2pitch(freq);

if frame2volume(frame)<ptOpt.volTh
	if plotOpt, fprintf('音量太小，音高設定為 0！\n'); end
	pitch=0;
end
if localMinCount>=ptOpt.maxAmdfLocalMinCount
	if plotOpt, fprintf('Local min. count = %d >= %d，音高設定為 0！\n', localMinCount, ptOpt.maxAmdfLocalMinCount); end
	pitch=0;
end

% ====== Plot related information
if plotOpt,
	clf;
	plotNum=3;
	frameAxisH=subplot(plotNum,1,1);
	frameH=plot(1:length(frame), frame, '.-'); axis([1, length(frame), -2^ptOpt.nbits/2, 2^ptOpt.nbits/2]); grid on; title('Frame');
	acfAxisH=subplot(plotNum,1,2);
	acfH=plot(1:length(acf), acf, '.-'); axis tight; title('ACF vector');
	amdfAxisH=subplot(plotNum,1,3);
	amdfH=plot(1:length(amdf), amdf, '.-'); axis tight; title('AMDF vector');
	if ~isempty(localMinIndex)
		localMinIndexH=line(localMinIndex, amdf(localMinIndex), 'linestyle', 'none', 'color', 'k', 'marker', 'o');
	end
	if ~isempty(minIndex)
		minIndexH=line(minIndex, amdf(minIndex), 'linestyle', 'none', 'color', 'r', 'marker', 'o');
		manualBarH=line(minIndex*[1 1], get(amdfAxisH, 'ylim'), 'color', 'm');
	end
	line(beginIndex*[1 1], get(gca, 'ylim'), 'color', 'r');
	line(  endIndex*[1 1], get(gca, 'ylim'), 'color', 'r');
	
	userDataSet;
	% 設定滑鼠按鈕的反應動作，以便讓使用者修正 pitch
	if gcf~=1;
		frameWinMouseAction;
	end
end

% ====== selfdemo
function selfdemo
waveFile='倫敦鐵橋垮下來_不詳_0.wav';
[y, fs, nbits]=wavReadInt(waveFile);
frameMat=buffer2(y, 256, 0);
frame=frameMat(:, 220);
ptOpt=ptOptSet(fs, nbits);
volume=frame2volume(frameMat);
ptOpt.volTh=getVolTh(volume, ptOpt.frameSize, 4);
plotOpt=1;
feval(mfilename, frame, plotOpt, ptOpt);

figure;
waveFile='宣王憤起揮天戈.wav';
[y, fs, nbits]=wavReadInt(waveFile);
frameMat=buffer2(y, 512, 0);
frame=frameMat(:, 33);
ptOpt=ptOptSet(fs, nbits);
ptOpt.frameSize=512;
ptOpt.overlap=0;
volume=frame2volume(frameMat);
ptOpt.volTh=getVolTh(volume, ptOpt.frameSize, 4);
plotOpt=1;
feval(mfilename, frame, plotOpt, ptOpt);