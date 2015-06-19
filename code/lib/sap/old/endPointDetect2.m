function [endPoint, voicedSegment] = endPointDetect(y, fs, nbits, plotOpt, epdPrm)
% endPointDetect: End point detection based on volume and zero-crossing rate
%	Usage: [endPoint, voicedSegment] = endPointDetect(y, fs, nbits, plotOpt, epdPrm)

%	Roger Jang, 20041118

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, nbits=8; end
if nargin<4, plotOpt=0; end
if nargin<5,
	epdPrm.frameSize = round(fs/31.25);	% fs=8000 ===> frameSize=256;
	epdPrm.frameSize = 320;		% same as htk
	epdPrm.overlap = round(epdPrm.frameSize*3/4);
	epdPrm.zcrRatio = 0.3;
	epdPrm.testFrameNum = 5;		% 測試平均環境噪音的音框數（沒用到）
	epdPrm.volRatio1=2;		
	epdPrm.volRatio2=4;
end

frameSize=epdPrm.frameSize;
overlap=epdPrm.overlap;
zcrRatio=epdPrm.zcrRatio;

% ====== Zero adjusted
y = y-round(mean(y));

% ====== Take frames
framedY  = buffer2(y, frameSize, overlap);
frameNum = size(framedY, 2);			% Number of frames

% ====== Compute volume and zcr
volume=zeros(frameNum, 1);
zcr=zeros(frameNum, 1);
for i=1:frameNum
	frame=framedY(:, i);
	frame=frame-round(mean(frame));
	volume(i)=sum(abs(frame-mean(frame)))/frameNum;
	zcr(i)=sum(abs(diff(frame>0)));
end
frameTime = frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
volTh=max(volume);
volTh1=epdPrm.volRatio1*volTh;
volTh2=epdPrm.volRatio2*volTh;

% ====== Find initial end points according volume level2 (upper level)
voicedIndex = find(volume>=volTh2);
while isempty(voicedIndex)	% 若找不到有聲音部分，改變音量門檻值
	volTh1=max(volume)/16;
	volTh2=max(volume)/8;
	voicedIndex = find(volume>=volTh2);
end

zcrTh = max(zcr)*zcrRatio;
sound = [];
k = 1;
sound(k).begin = voicedIndex(1);
for i=2:length(voicedIndex)-1,
	if voicedIndex(i+1)-voicedIndex(i)>1,
		sound(k).end = voicedIndex(i);
		sound(k+1).begin = voicedIndex(i+1);
		k = k+1;
	end
end
sound(k).end = voicedIndex(end);

% ====== Delete short sound clips
%index = [];
%for i=1:length(sound),
%	if (sound(i).end-sound(i).begin)<4
%		index = [index, i];
%	end
%end
%sound(index) = [];

% ====== Expand end points to volume level1 (lower level)
for i=1:length(sound),
	head = sound(i).begin;
	while (head-1)>=1 & volume(head-1)>volTh1,
		head=head-1;
	end
	sound(i).begin = head;
	tail = sound(i).end;

	while (tail+1)<=length(volume) & volume(tail+1)>volTh1,
		tail=tail+1;
	end
	sound(i).end = tail;
end

% ====== Expand end points to include high zcr region
for i=1:length(sound),
	head = sound(i).begin;
	while (head-1)>=1 & zcr(head-1)>zcrTh,
		head=head-1;
	end
	sound(i).begin = head;
end

% ====== Delete repeated sound segments
if length(sound) ~=0,
	index = [];
	for i=1:length(sound)-1,
		if sound(i).begin==sound(i+1).begin & sound(i).end==sound(i+1).end,
			index=[index, i];
		end
	end
	sound(index) = [];
end;

% ====== Transform sample-point-based index
if length(sound) ~=0,
	for i=1:length(sound),
		voicedSegment(i).begin = sound(i).begin*(frameSize-overlap);
		voicedSegment(i).end   = min(length(y), sound(i).end*(frameSize-overlap));
	end
	endPoint=[voicedSegment(1).begin, voicedSegment(end).end];	% 取頭尾
else
	endPoint = [];
end;

if plotOpt,
	axes1H=subplot(4,1,1);
	time=(1:length(y))/fs;
	plot(time, y); grid on
	axis([min(time), max(time), -2^nbits/2, 2^nbits/2]);
	ylabel('Amplitude');
	title('Waveform');
	% Plot end points
	yBound=[-2^nbits/2, 2^nbits/2];
	for i=1:length(sound),
		line(frame2sampleIndex(sound(i).begin, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  sound(i).end, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes2H=subplot(4,1,2);
	plot(frameTime, volume, '.-'); grid on
	line([min(frameTime), max(frameTime)], volTh1*[1 1], 'color', 'c');
	line([min(frameTime), max(frameTime)], volTh2*[1 1], 'color', 'c');
	axis tight
	ylabel('Volume');
	title('Volume');
	% Plot end points
	yBound = [min(volume) max(volume)];
	for i=1:length(sound),
		line(frame2sampleIndex(sound(i).begin, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  sound(i).end, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes3H=subplot(4,1,3);
	plot(frameTime, zcr, '.-'); grid on
	line([min(frameTime), max(frameTime)], zcrTh*[1 1], 'color', 'c');
	axis([min(frameTime), max(frameTime), 0, max(zcr)]);
	ylabel('ZCR');
	title('Zero crossing rate');
	% Plot end points
	yBound = [0 max(zcr)];
	for i=1:length(sound),
		line(frame2sampleIndex(sound(i).begin, frameSize, overlap)/fs*[1,1], yBound, 'color', 'm');
		line(frame2sampleIndex(  sound(i).end, frameSize, overlap)/fs*[1,1], yBound, 'color', 'g');
	end

	axes4H=subplot(4,1,4);
	voicedIndex=endPoint(1):endPoint(2);
	voicedTime=time(voicedIndex);
	voicedY=y(voicedIndex);
	voicedH=plot(voicedTime, voicedY); grid on
	axis([time(endPoint(1)), time(endPoint(2)), -inf, inf]);
	ylabel('Amplitude');
	title('Voiced waveform');
	
	U.y=y; U.fs=fs; U.nbits=nbits;
	U.axes1H=axes1H; U.axes2H=axes2H; U.axes3H=axes3H; U.axes4H=axes4H;
	U.voicedIndex=voicedIndex; U.voicedH=voicedH;
	U.voicedY=voicedY; U.voicedTime=voicedTime;
	set(gcf, 'userData', U);
	
	uicontrol('string', 'Play all', 'callback', 'U=get(gcf, ''userData''); sound(U.y/(2^U.nbits/2), U.fs);');
	uicontrol('string', 'Play voiced', 'callback', 'U=get(gcf, ''userData''); sound(U.voicedY/(2^U.nbits/2), U.fs);', 'position', [100, 20, 60, 20]);

	% Play the segmented sound
%	head = sound(1).begin*(frameSize-overlap);
%	tail = min(length(y), sound(end).end*(frameSize-overlap));
%	thisY = y(head:tail);
%	fprintf('His return to hear the cutted sound %g:', i);
%	pause;
%	fprintf('\n');
%	wavplay(thisY, fs, 'sync');
%	fprintf('\n');
end

% ======
function sampleIndex=frame2sampleIndex(frameIndex, frameSize, overlap)
sampleIndex=(frameIndex-1)*(frameSize-overlap)+round(frameSize/2);
	
% ====== Self demo
function selfdemo
waveFile='malisa\SenPC0000_2.wav';
waveFile='Arthas\SenIC0000_3.wav';
waveFile='tracy\SenIC0000_2.wav';
waveFile='abo\SenIC0000_3.wav';
waveFile='roger\SenIC0000_2.wav';
waveFile='jacky\SenIC0000_3.wav';
waveFile='__Cyberon__\2_f_0.wav';
waveFile='__Cyberon__\1_f_19.wav';
waveFile='主人下馬客在船.wav';
waveFile='此恨綿綿無絕期+sil.wav';
waveFile='楊家有女初長成.wav';
plotOpt = 1;
[y, fs, nbits] = wavReadInt(waveFile);
endPoint = feval(mfilename, y, fs, nbits, plotOpt);