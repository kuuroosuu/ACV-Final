function candidate=amdf2pitchCandidate(amdf, plotOpt)

if nargin<1, selfdemo; return; end
if nargin<2, plotOpt=0; end

fs=8000;
origAmdf=amdf;
maxValue=max(amdf);
%amdf(1:7)=maxValue;	% 前 7 點不予考慮 local minima
%minValue=min(amdf);
%range=maxValue-minValue;

localMinIndex=find(localMax(-amdf));
localMinIndex(1)=[];	% 第一點不考慮
if localMinIndex(end)==length(amdf)
	localMinIndex(end)=[];	% 最後一點如果是 amdf 的最後一點，也不考慮
end
localMinValue=amdf(localMinIndex);
% 只留 n 點最小點
n=10;
[junk, index]=sort(localMinValue);
if length(index)>n
	threshold=localMinValue(index(n));
	index=find(localMinValue>threshold);
	localMinValue(index)=[];
	localMinIndex(index)=[];
end

prob=(maxValue-localMinValue)/maxValue;
prob=prob/sum(prob);
prob=log(prob);

for i=1:length(prob)
	candidate(i).index=localMinIndex(i);
	candidate(i).stateProb=prob(i);
	freq=10*floor((fs+floor((localMinIndex(i)-1)/2))/(localMinIndex(i)-1));
	candidate(i).pitch=freq2pitch(freq);
%	freq=fs/(localMinIndex(i)-1);
%	candidate(i).pitch=10*freq2semitone(freq);
end

if plotOpt
%	subplot(2,1,1);
%	plot(origAmdf, '.-');
%	set(gca, 'xlim', [-inf inf]);
%	title('Original AMDF');
%	subplot(2,1,1);
	plot(amdf, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('AMDF vector for candidate selection (top figure: pitch, bottom figure: log prob.');
	line(localMinIndex, localMinValue, 'color', 'r', 'marker', 'o', 'linestyle', 'none');
	grid on;
	for i=1:length(candidate),
		h=text(candidate(i).index, amdf(candidate(i).index), num2str(candidate(i).pitch));
		set(h, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
		h=text(candidate(i).index, amdf(candidate(i).index), num2str(candidate(i).stateProb));
		set(h, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
	end
end

% ====== Self demo
function selfdemo
waveFile='../spce4singing/waveData/cweb/frida_faraway/this old man.wav';
%waveFile='waveData/cweb/roger_faraway/the more we get together.wav';
%waveFile='../spce4singing/waveData/cweb/beball/the more we get together.wav';
[y, fs, nbits]=waveFileRead(waveFile);
[PP, CP]=setParam;
[pitch, volume, PP] = wave2pitchVolume(y, 15, 0, PP);
volTh=getVolumeThreshold(volume);
index=find(volume<volTh);
framedY=buffer2(y, PP.frameSize, PP.overlap);
framedY(:, index)=[];
frameIndex=18;
frame=framedY(:, frameIndex);
frame2=frameFlip(frame);
frame3=localAverage(frame2);
amdf=frame2amdf(frame3, 128);
plotOpt=1;
feval(mfilename, amdf, plotOpt);
xlabel(sprintf('Frame %d of "%s"', frameIndex, waveFile));