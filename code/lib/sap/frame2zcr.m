function zcr = frame2zcr(frameMat, method, shiftAmount);
% frame2zcr: Frame or frame matrix to ZCR (zero-crossing rate) conversion
%	Usage: zcr = frame2zcr(frame, method, shiftAmount)
%			frame is a column vector
%	       zcr = frame2zcr(frameMat, method, shiftAmount);
%			frameMat is a matrix with columns of frame
%
%		method: 1 for "zero-value samples are not counted for zcr"
%			2 for "zero-value samples are counted for zcr"
%
%		Note the framewise zero-justification is performed for each frame.
%
%		Regarding the selection of shiftAmount, please refer to the following example.
%
%	For example:
%		waveFile='csNthu.wav';
%		[y, fs, nbits]=wavreadInt(waveFile);
%		frameSize=320; overlap=160;
%		framedY=buffer2(y, frameSize, overlap);
%		frameNum=size(framedY, 2);
%		volume=frame2volume(framedY);
%		[minVolume, index]=min(volume);
%		shiftAmount=2*max(abs(framedY(:,index)));	% twice the max amplitude of the minimum-volume frame
%		zcr1=frame2zcr(framedY, 1, 0);
%		zcr2=frame2zcr(framedY, 1, shiftAmount);
%		time=(1:length(y))/fs;
%		frameTime=((0:frameNum-1)*(frameSize-overlap)+0.5*frameSize)/fs;
%		subplot(2,1,1); plot(time, y); axis([min(time), max(time), 2^nbits/2*[-1 1]]); ylabel('Amplitude'); title(waveFile);
%		subplot(2,1,2); plot(frameTime', [zcr1; zcr2]', '.-'); axis tight; ylabel('ZCR'); legend('shiftAmount=0', ['shiftAmount=', num2str(shiftAmount)]);

%	Roger Jang, 20060303, 20070417

if nargin<1, selfdemo; return; end
if nargin<2, method=1; end
if nargin<3, shiftAmount=0; end

[frameSize, frameNum]=size(frameMat);

% Zero justification
for i=1:frameNum
%	frameMat(:,i)=frameMat(:,i)-round(mean(frameMat(:,i)));
	frameMat(:,i)=frameMat(:,i)-mean(frameMat(:,i));
end

zcr=zeros(1, frameNum);
switch method
	case 1
		for i=1:frameNum
			frame=frameMat(:, i)-shiftAmount;
			zcr(i)=sum(frame(1:end-1).*frame(2:end)<0);	% signals at zero do not count for zcr
		end
	case 2
		for i=1:frameNum
			frame=frameMat(:, i)-shiftAmount;
			zcr(i)=sum(frame(1:end-1).*frame(2:end)<=0);	% signals at zero do count for zcr
		end
	otherwise
		error('Unknown method!');
end

% ====== Self demo
function selfdemo
waveFile='csNthu8k.wav';
[y, fs, nbits]=wavread(waveFile);
frameSize=256; overlap=0;
frameMat=buffer2(y, frameSize, overlap);
frameNum=size(frameMat, 2);
zcr1=frame2zcr(frameMat, 1);
zcr2=frame2zcr(frameMat, 2);
time=(1:length(y))/fs;
frameTime=((0:frameNum-1)*(frameSize-overlap)+0.5*frameSize)/fs;
%subplot(2,1,1); plot(time, y); axis([min(time), max(time), 2^nbits/2*[-1 1]]); ylabel('Amplitude'); title(waveFile);
subplot(2,1,1); plot(time, y); axis([min(time), max(time), [-1 1]]); ylabel('Amplitude'); title(waveFile);
subplot(2,1,2); plot(frameTime', [zcr1; zcr2]', '.-'); axis tight; ylabel('ZCR'); legend('Method 1', 'Method 2');

figure
waveFile='csNthu.wav';
[y, fs, nbits]=wavread(waveFile);
frameSize=320; overlap=160;
framedY=buffer2(y, frameSize, overlap);
frameNum=size(framedY, 2);
volume=frame2volume(framedY);
[minVolume, index]=min(volume);
shiftAmount=2*max(abs(framedY(:,index)));	% twice the max amplitude of the minimum-volume frame
zcr1=frame2zcr(framedY, 1, 0);
zcr2=frame2zcr(framedY, 1, shiftAmount);
time=(1:length(y))/fs;
frameTime=((0:frameNum-1)*(frameSize-overlap)+0.5*frameSize)/fs;
%subplot(2,1,1); plot(time, y); axis([min(time), max(time), 2^nbits/2*[-1 1]]); ylabel('Amplitude'); title(waveFile);
subplot(2,1,1); plot(time, y); axis([min(time), max(time), [-1 1]]); ylabel('Amplitude'); title(waveFile);
subplot(2,1,2); plot(frameTime', [zcr1; zcr2]', '.-'); axis tight; ylabel('ZCR'); legend('shiftAmount=0', ['shiftAmount=', num2str(shiftAmount)]);