addpath d:/users/jang/matlab/toolbox/utility
addpath d:/users/jang/matlab/toolbox/audioProcessing

waveFile='bug_1.wav';
[y, fs, nbits]=wavread(waveFile);
len=length(y)/fs;
fprintf('Length of wave = %f sec\n', len);
y1=y(:, 1);
y2=y(:, 2);
%subplot(2,1,1); plot(y1);
%subplot(2,1,2); plot(y2);
tic
minDuration=5;
plotOpt=0;
cutPointSampleIndex=waveSegment(y2, fs, minDuration, plotOpt);
fprintf('Time for waveSegment = %f sec\n', toc);
fprintf('No. of segments = %d\n', length(cutPointSampleIndex)+1);

% Playback of segments
temp=[1, cutPointSampleIndex, length(y)];
segmentNum=length(temp)-1;
segmentDuration=diff(temp)/fs;
fprintf('Duration of segments = %s\n', mat2str(segmentDuration));
for i=1:segmentNum
	thisWave=y2(temp(i):temp(i+1));
	fprintf('Press any key to play segment %d/%d (%f sec)...', i, segmentNum, segmentDuration(i)); pause; fprintf('\n');
	wavplay(thisWave, fs, 'sync');
end