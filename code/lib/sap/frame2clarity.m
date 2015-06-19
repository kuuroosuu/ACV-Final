function [clarity, pfMat] = frame2clarity(frameMat, fs, pfType, pfMethod, plotOpt);
% frame2clarity: Frame (or frame matrix) to clarity
%
%	Usage:
%		[clarity, pf] = frame2clarity(frameMat, fs, pfType, pfMethod, plotOpt);
%			frameMat: a column vector of a frame, or a matrix where each column is a frame
%			fs: sampling rate
%			pfType: 'acf', 'nsdf', or 'amdf'
%
%	Example:
%		waveFile='二十四橋明月夜.wav';
%		[y, fs, nbits]=wavread(waveFile);
%		frameSize=512; overlap=0;
%		frameMat=enframe(y, frameSize, overlap);
%		frameMat=frameZeroMean(frameMat, 2);
%		frameNum=size(frameMat, 2);
%		clarity1=frame2clarity(frameMat, fs, 'acf');
%		clarity2=frame2clarity(frameMat, fs, 'nsdf');
%		clarity3=frame2clarity(frameMat, fs, 'amdf');
%		time=(1:length(y))/fs;
%		frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
%		subplot(2,1,1); plot(time, y); axis tight; ylabel('Amplitude'); title(waveFile);
%		subplot(2,1,2); plot(frameTime, [clarity1; clarity2; clarity3]', '.-'); axis tight; ylabel('Clarity'); grid on
%		legend('pfType=acf', 'pfType=nsdf', 'pfType=amdf');
%		xlabel('Time (sec)');

%	Roger Jang, 20070417

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, pfType='acf'; end
if nargin<4, pfMethod=2; end
if nargin<5, plotOpt=0; end

[frameSize, frameNum]=size(frameMat);
% Zero justification
for i=1:frameNum
	frameMat(:,i)=frameMat(:,i)-mean(frameMat(:,i));
end

minPitch=30;
maxPitch=84;
index1=round(fs/pitch2freq(maxPitch));
index2=min(round(fs/pitch2freq(minPitch)), frameSize);

clarity=zeros(1, frameNum);
pfMat=zeros(frameSize, frameNum);
switch pfType
	case 'acf'
		for i=1:frameNum
			frame=frameMat(:,i);
			pf=frame2acfMex(frame, length(frame), pfMethod);
			pfMat(:,i)=pf;
			maxPf=max(pf(index1:index2));
			if pf(1)==0
				clarity(i)=0;
			else
				clarity(i)=maxPf/pf(1);
			end
		end
	case 'nsdf'
		for i=1:frameNum
			frame=frameMat(:,i);
			pf=frame2nsdfMex(frame, length(frame), pfMethod);
			pfMat(:,i)=pf;
			clarity(i)=max(pf(index1:index2));
		end
	case 'amdf'
		for i=1:frameNum
			frame=frameMat(:,i);
			pf=frame2amdfMex(frame, length(frame), pfMethod);
			pfMat(:,i)=pf;
			minAmdf=min(pf(index1:index2));
			maxAmdf=max(pf(index1:index2));
			if maxAmdf==0
				clarity(i)=0;
			else
				clarity(i)=1-minAmdf/maxAmdf;
			end
		end
	otherwise
		error('Unknown pfType!');
end

if plotOpt
	if frameNum==1
		subplot(2,1,1); plot(frame); set(gca, 'xlim', [-inf inf]);
		subplot(2,1,2); plot(pf); set(gca, 'xlim', [-inf inf]);
		axisLimit=axis;
		line(index1*[1 1], axisLimit(3:4), 'color', 'r');
		line(index2*[1,1], axisLimit(3:4), 'color', 'r');
	end
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);