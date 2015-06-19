function flag=waveVolumeCheck(wObj, volCheckPrm, plotOpt)
% waveVolumeCheck: Check the volume of a wave object
%	Usage: flag=waveVolumeCheck(wObj)
%
%		flag=0 ==> 正常
%		flag=1 ==> 音量太小
%		flag=2 ==> 環境雜訊太大
%		flag=3 ==> 爆音太多

%	Roger Jang, 20101213

if nargin<1, selfdemo; return; end

if ischar(wObj)
	wObj=waveFile2obj(wObj);
end

if nargin<2 | isempty(volCheckPrm)
	volCheckPrm.maxVolTh=140;
	volCheckPrm.minVolTh=50;
	volCheckPrm.satAmpRatio=0.95;
	volCheckPrm.satSampleCountRatioTh=0.005;
end
if nargin<3, plotOpt=0; end

if wObj.amplitudeNormalized
	wObj.signal=wObj.signal*2^wObj.nbits/2;
end
frameSize=320;
overlap=160;
largeCount=0;
smallCount=0;
method=1;
flag=0;			% 預設正常，不需拒絕

volume=wave2volume(wObj.signal, frameSize, overlap, method);
frameNum=length(volume);
maxVol=max(volume);
minVol=min(volume);
if (maxVol<volCheckPrm.maxVolTh*frameSize)	% 音量太小
	flag=1;
end

if (minVol>volCheckPrm.minVolTh*frameSize)	% 使用全向式麥克風，或是環境雜訊太大
	flag=2;
end

ampTh=(2^wObj.nbits/2)*volCheckPrm.satAmpRatio;
largeCount=sum(wObj.signal> ampTh);
smallCount=sum(wObj.signal<-ampTh);
if (largeCount+smallCount)/length(volume)>volCheckPrm.satSampleCountRatioTh
	flag=3;
end

if plotOpt
	subplot(2,1,1);
	time=(1:length(wObj.signal))/wObj.fs;
	plot(time, wObj.signal);
	line([min(time), max(time)],  ampTh*[1 1], 'color', 'r');
	line([min(time), max(time)], -ampTh*[1 1], 'color', 'r');
	axis([-inf inf 2^wObj.nbits/2*[-1 1]]);
	subplot(2,1,2);
	plot(1:frameNum, volume);
	line([1, frameNum], frameSize*volCheckPrm.minVolTh*[1 1], 'color', 'r');
	line([1, frameNum], frameSize*volCheckPrm.maxVolTh*[1 1], 'color', 'g');
	maxVolLimit=max(1000000, maxVol);
	axis([-inf inf -inf maxVolLimit]);
	legend('Volume', 'minVolTh', 'maxVolTh');
end

% ====== Self demo
function selfdemo
waveFile='yi_cuen_xiang_s_yi_cuen_huei.wav';
flag=waveVolumeCheck(waveFile, [], 1);