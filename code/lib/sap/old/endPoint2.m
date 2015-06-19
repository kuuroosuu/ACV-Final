function output=endPoint(y, fs, nbits, epdPrm, plotOpt)
% endPoint: Endpoint detection（適用於 RMC 的電話語音錄音）
%	Usage: 	output=endPoint(y, fs, nbits, epdPrm, plotOpt)
%		output=endPoint(waveFile, [], [], epdPrm, plotOpt)
%
%		output(1): start index
%		output(2): end index

if nargin==0; selfdemo; return; end
if isstr(y); file=y; [y, fs, nbits]=wavread(file); end
if nargin<4 | isempty(epdPrm)
	epdPrm.duDuration=0.2;	% 開頭部分的「嘟」聲，約 0.2 秒
	epdPrm.volRatio=0.15;
	epdPrm.zcrRatio=0.35;
	epdPrm.frameSize=256;
	epdPrm.overlap=128;
	epdPrm.extension=2;		% 從原始端點向前後各延伸2個frame
end
if nargin<5; plotOpt=0; end

y=y*2^nbits/2;		% 轉成整數
%y=y((epdPrm.duDuration*fs+1):end);	% 刪除電話語音錄音「嘟」的部分
startFrameIndex=1+floor((epdPrm.duDuration*fs-epdPrm.overlap)/(epdPrm.frameSize-epdPrm.overlap));	% 刪除電話語音錄音「嘟」的部分
%y=y-mean(y);

time=(1:length(y))/fs;
framedY=buffer2(y, epdPrm.frameSize, epdPrm.overlap);
framedY=framedY-ones(epdPrm.frameSize,1)*mean(framedY);	% 設定每個 frame 的平均值為零
frameNum=size(framedY, 2);
%fprintf('frameNum=%d\n', frameNum);
frameTime=((0:frameNum-1)*(epdPrm.frameSize-epdPrm.overlap)+epdPrm.frameSize/2)/fs;

% ====== Compute volume
vol=sum(abs(framedY));
vol(1:startFrameIndex-1)=0;
%disp(vol(1:3));
volTh=max(vol)*epdPrm.volRatio;		% Volume threshold
%fprintf('volTh=%g\n', volTh);

% ====== Compute zero crossing rate
zcr=sum(abs(diff(framedY>0)));
zcr=sum(framedY(1:end-1,:).*framedY(2:end, :)<0);
zcr(1:startFrameIndex-1)=0;
%disp(zcr(1:3));
zcrTh=max(zcr)*epdPrm.zcrRatio;
%fprintf('zcrTh=%g\n', zcrTh);

% ====== Compute endpoint
index1=vol>volTh;	% 音量要夠大
index2=zcr<zcrTh;	% 過零率要夠小
index=index1 & index2;
% 找出從開始連續不為零的 index
firstNonZero=[];
for i=startFrameIndex:length(index)
	if index(i)>0
		firstNonZero=[firstNonZero, i];
	else
		break;
	end
end
index(firstNonZero)=0;

% 在第一秒內，找出第一個 [0 1 0]
for i=startFrameIndex:ceil(fs/(epdPrm.frameSize-epdPrm.overlap))
	if all(index(i:i+2)==[0 1 0])
		break;
	end
end
patternIndex1=i+1;
index(patternIndex1)=0;

% 在第一秒內，找出第一個 [0 1 1 0]
for i=startFrameIndex:ceil(fs/(epdPrm.frameSize-epdPrm.overlap))
	if all(index(i:i+3)==[0 1 1 0])
		break;
	end
end
patternIndex2=[i+1, i+2];
index(patternIndex2)=0;

temp=find(index);
%fprintf('start frame = %d\n', temp(1));
%fprintf('end frame = %d\n', temp(end));

startIndex=(temp(1)-1)*(epdPrm.frameSize-epdPrm.overlap)+epdPrm.frameSize/2-epdPrm.extension*epdPrm.frameSize;
endIndex=(temp(end)-1)*(epdPrm.frameSize-epdPrm.overlap)+epdPrm.frameSize/2+epdPrm.extension*epdPrm.frameSize;
output(1)=max(startIndex, 1);
output(2)=min(endIndex, length(y));
%fprintf('start sample = %d\n', output(1));
%fprintf('end sample = %d\n', output(2));

% ====== Code for plotting
if plotOpt
	plotNum=3;
	subplot(plotNum,1,1);
	plot(time, y);
	set(gca, 'xlim', [min(time), max(time)]); grid on
	if exist('file'), title(file); end

	subplot(plotNum,1,2);
	plot(frameTime, vol, '.-');
	line([min(frameTime), max(frameTime)], volTh*[1 1], 'color', 'r');
	set(gca, 'xlim', [min(frameTime), max(frameTime)]); grid on
	ylabel('Volume');

	subplot(plotNum,1,3);
	plot(frameTime, zcr, '.-');
	line([min(frameTime), max(frameTime)], zcrTh*[1 1], 'color', 'r');
	set(gca, 'xlim', [min(frameTime), max(frameTime)]); grid on
	ylabel('Zero crossing rate');

	subplot(plotNum, 1, 1);
	limit=axis;
	line(output(1)*[1 1]/fs, limit(3:4), 'color', 'r');
	line(output(2)*[1 1]/fs, limit(3:4), 'color', 'r');
	line(epdPrm.duDuration*[1 1], limit(3:4), 'color', 'm');
	subplot(plotNum, 1, 2);
	line(frameTime(index), vol(index), 'marker', 'o', 'color', 'r', 'linestyle', 'none');
	
	line(frameTime(firstNonZero), vol(firstNonZero), 'marker', 'x', 'color', 'k', 'linestyle', 'none');	% 被移除的點：開始連續非零
	line(frameTime(patternIndex1), vol(patternIndex1), 'marker', 'x', 'color', 'k', 'linestyle', 'none');	% [0 1 0]
	line(frameTime(patternIndex2), vol(patternIndex2), 'marker', 'x', 'color', 'k', 'linestyle', 'none');	% [0 1 1 0]
	
	newY=y(output(1):output(2));
	sound(newY/(2^nbits/2), fs);
return

Y=fft(framedY);
z=log(Y.*conj(Y)+eps);
z=Y.*conj(Y);
z=z(1:65, :);
%z=z*diag(1./sum(z));	% Sum of each column is 1
variance=diag(cov(z));
th=max(variance)/5;
subplot(plotNum,1,4);
plot(frameTime, variance, '.-');
set(gca, 'xlim', [min(frameTime), max(frameTime)]); grid on
line([min(frameTime), max(frameTime)], th*[1 1], 'color', 'r');
%set(gca, 'ylim', [0 20000]);
ylabel('Spectral variance');

temp=abs(framedY)+eps;
temp=temp*diag(1./sum(temp));	% Sum of each column is 1
entropy=sum(temp.*log(temp));
subplot(plotNum,1,5);
plot(frameTime, entropy, '.-');
set(gca, 'xlim', [min(frameTime), max(frameTime)]); grid on
ylabel('Spectral entropy');

diffVal=sum(abs(diff(framedY)));
subplot(plotNum,1,6);
plot(frameTime, diffVal, '.-');
set(gca, 'xlim', [min(frameTime), max(frameTime)]); grid on
ylabel('diffVal');
end

% ====== self demo
function selfdemo
waveFile='alexxx1.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\car3.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\coco5.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\eva2.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\irene3.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\coco1.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\coco4.wav';
waveFile='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\soph5.wav';
plotOpt=1;
feval(mfilename, waveFile, [], [], [], plotOpt);