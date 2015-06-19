function [pitch, volume, volTh, amdfMat, A]=wave2pitchByAmdfDp(y, fs, nbits, PP, plotOpt)
% wave2pitchByDp: 利用 DP 的方法來找出最佳的音高曲線
%	Usage: [pitch, amdfMat, A]=wave2pitchByDp(y, fs, nbits, PP, plotOpt)
%
%	這是一個簡化的版本，假設送進來的y都是有聲音的母音（都有音高），效果會比較好
%	如果送進來的 wave 其中有任一個音框沒有音高，amdf 可能不出現 local min.，本函數將會分段進行 DP。
%	此外，此函數的結果會根據 PP.indexDiffWeight 而有不同的結果，最好先經大量測試，取得此參數的最佳值。

%	Roger Jang, 20050228, 20051019

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, nbits=8; end
if nargin<4, PP=ptOptSet(fs, nbits); end
if nargin<5, plotOpt=0; end

frameMat=buffer2(y, PP.frameSize, PP.overlap);
volume=frame2volume(frameMat);
volTh=max(volume)/10;
amdfMat=frameMat2amdfMat(frameMat, PP.frameSize/2, 3);
[amdfLen, frameNum]=size(amdfMat);
if frameNum==1, error('Error: Need at least two frames!'); end

% 產生 A 結構，記錄每一個音框的 PPC
for i=1:frameNum
	% ====== 計算每一個音框的 PPC (Pitch Period Candidate)
	[ppcIndex, ppcValue] = amdf2ppc(amdfMat(:,i), fs, nbits, PP);		% ppcValue 和 amdf(ppcIndex) 可能不同，因為可能經過倍頻轉換的修改
	% ====== 儲存資料於 A 結構
	A(i).ppcLen=length(ppcIndex);
	for j = 1:length(ppcIndex)
		A(i).ppc(j).index = ppcIndex(j);
		A(i).ppc(j).amdfValue = ppcValue(j);
		A(i).ppc(j).totalDist=0;
		A(i).ppc(j).from=-1;
		A(i).pp=0;
	end
end

% 先進行分段，再對每段使用 DP 進行音高追蹤
headPos = 1;
pitch = zeros(1, frameNum);
while headPos < frameNum
	while A(headPos).ppcLen == 0
		headPos = headPos + 1;
		if headPos == frameNum
			break;
		end
	end
	
	tailPos = headPos;
	while A(tailPos).ppcLen ~= 0
		tailPos = tailPos + 1;
		if tailPos == frameNum + 1
			break;
		end
	end
	
	range=headPos:tailPos-1;
	if length(range)>0
		[pitch(range), A(range)]=A2pitchByDp(A(range), fs, PP);
	end
	
	headPos = tailPos;
end

% 將 A 寫入檔案
%fid=fopen('a2.txt', 'w');
%for i=1:frameNum
%	fprintf(fid, 'A[%d].ppcLen=%d\n', i, A(i).ppcLen);
%	for j=1:A(i).ppcLen
%		fprintf(fid, '\tppc[%d]: amdfValue=%d, totalDist=%d, from=%d, index=%d\n', j, A(i).ppc(j).amdfValue, A(i).ppc(j).totalDist, A(i).ppc(j).from, A(i).ppc(j).index);
%	end	
%end
%fclose(fid);

pitch2=pitch;
if PP.useEpd
	pitch(volume<volTh)=0;		% 只保留高音量的部分
end

if plotOpt
	plotNum=4;
	frameTime=frame2sampleIndex(1:frameNum, PP.frameSize, PP.overlap)/fs;
	% ====== Plot wave
	waveAxisH=subplot(plotNum,1,1);
	time=(1:length(y))/fs;
	waveH=plot(time, y-mean(y));
	set(gca, 'xlim', [min(time), max(time)]);
	ylabel('Waveform');
	lFrameH=line(nan*[1 1], get(waveAxisH, 'ylim'), 'erase', 'xor', 'color', 'r');
	rFrameH=line(nan*[1 1], get(waveAxisH, 'ylim'), 'erase', 'xor', 'color', 'r');
	% ====== Plot volume
	subplot(plotNum,1,2);
	plot(frameTime, volume, '.-'); ylabel('Volume'); grid on;
	set(gca, 'xlim', [min(frameTime), max(frameTime)]); set(gca, 'ylim', [-inf inf]);
	line([min(frameTime), max(frameTime)], volTh*[1 1], 'color', 'r');
	% ====== Plot AMDF map
	subplot(plotNum,1,3);
	amdfPlot(amdfMat, frameTime, A); ylabel('AMDF matrix');
	barH(3)=line(nan*[1 1], get(gca, 'ylim'), 'erase', 'xor', 'color', 'r');
	% ====== Plot the final pitch
	subplot(plotNum,1,4);
	temp=pitch; temp(temp==0)=nan;
	temp2=pitch2; temp2(temp2==0)=nan;
	plot(frameTime, temp2, '.-', frameTime, temp, 'ro'); ylabel('Final pitch'); grid on;
	set(gca, 'xlim', [min(frameTime), max(frameTime)]); set(gca, 'ylim', [-inf inf]);
	barH(5)=line(nan*[1 1], get(gca, 'ylim'), 'erase', 'xor', 'color', 'r');
	% ====== Buttons for playback
	waveObj.signal=y-mean(y); waveObj.fs=fs; waveObj.nbits=nbits;
	pitchObj.signal=pitch; pitchObj.frameRate=fs/(PP.frameSize-PP.overlap);
	buttonH=wavePitchPlayButton(waveObj, pitchObj);
	% ====== For displaying frame2ppc
%	userDataSet;
%	set(gcf, 'WindowButtonDownFcn', 'DpWindowButtonDownFcn');
end

% ====== Selfdemo
function selfdemo
waveFile='lately2.wav';		% 16 kHz, 16 bits
waveFile='但使龍城飛將在.wav';
[y, fs, nbits]=wavReadInt(waveFile);
PP=ptOptSet(fs, nbits);
PP.useEpd=1;
PP.indexDiffWeight=10000;
%frameRange=[97, 134];		% Voiced range
%sampleRange=frame2sampleRange(frameRange, 512, 512-160);
%y=y(sampleRange(1):sampleRange(2));	% Voiced part
[pitch, amdfMat, A]=feval(mfilename, y, fs, nbits, PP, 1);
%save pitch pitch PP