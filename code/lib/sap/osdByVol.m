function timeSequence = osdByVol(y, fs, itera, plotOpt)
%osdByVol: Onset detection by volume
%	Usage: timeSequence = osdByVol(wave, fs, itera)

%	Hong-Ru Lee 20080820, Roger Jang 20080907

if nargin<1, selfdemo; return; end
if nargin<2, fs=11025; end
if nargin<3, itera=0.68; end
if nargin<4, plotOpt=1; end

wave=y;
wave=wave-mean(wave);

frameSize = 256;
overlap = 250;

% ====== high-pass filter
cutOffFreq=500;		% Cutoff freq
filterOrder=2;		% Order of filter 
[b, a]=butter(filterOrder, cutOffFreq/(fs/2), 'high');
wave=filter(b, a, y);

% ====== Set up parameters and find energy profile
frameMat = buffer2(wave, frameSize, overlap);
frameNum=size(frameMat, 2);
energy=sum(abs(frameMat));

% ====== Find local maxima
index1 = find(localMax(energy));
[sortedEnergy, inx]=sort(energy(index1));
[tmp, tmp1]=sort(diff(sortedEnergy,2));		% Find the biggest gap in energy
a1=min(sortedEnergy(tmp1(end-10:end)));
b1=sortedEnergy(tmp1(end));
if b1<a1
    tmpab = b1;
    b1 = a1;
    a1 = tmpab;
end
threshold = a1+itera*(b1-a1);
index2 = find(energy>threshold);

index = intersect(index1, index2);

% ====== Eliminate close candidates
index(find(diff(index)<80))=[];

timeSequence = frame2sampleIndex(index, frameSize, overlap)/fs;

if plotOpt==1
	subplot(3,1,1);
	time=(1:length(y))/fs;
	plot(time, y); axis tight
	subplot(3,1,2);
	plot(time, wave); axis tight
	subplot(3,1,3);
	frameTime = frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	plot(frameTime, energy); axis([min(frameTime), max(frameTime), min(energy), max(energy)])
	line(frameTime(index1), energy(index1), 'marker', '.', 'color', 'r', 'linestyle', 'none');
	line([min(frameTime), max(frameTime)], threshold*[1 1], 'color', 'k');
	for i=1:length(index)
		line(frameTime(index(i))*[1 1], [min(energy), max(energy)], 'color', 'm');
	end
	subplot(3,1,1);
	for i=1:length(timeSequence)
		line(timeSequence(i)*[1 1], [min(y), max(y)], 'color', 'm');
	end
end

% ====== Self demo
function selfdemo
waveFile='0385_noisyTapping.wav';
[y, fs, nbits]=wavReadInt(waveFile);
itera=0.5;
plotOpt=1;
timeSequence=osdByVol(y, fs, itera, plotOpt);