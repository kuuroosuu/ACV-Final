clear all; addMyPath;
fprintf('Compiling ptByDpOverPfMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility ptByDpOverPfMex.cpp /users/jang/c/lib/audio/audio.cpp /users/jang/c/lib/utility/utility.cpp
fprintf('Finish compiling ptByDpOverPfMex.cpp\n');
dos('copy /y ptByDpOverPfMex.mex* ..');
addpath ..
fprintf('Test ptByDpOverPfMex...\n');

waveFile='../twinkle_twinkle_little_star.wav';
waveFile='/dataSet/childSong/2007-音訊處理與辨識/19461108任佳王民/十個印第安人_不詳_0.wav';
pvFile=  [waveFile(1:end-4), '.pv'];
pv=asciiRead(pvFile); 

fprintf('waveFile=%s\n', waveFile);
[y, fs, nbits]=wavReadInt(waveFile);
pfType=1; ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.useVolThreshold=1;
ptOpt.cutLeadingTrailingZero=0;
ptOpt.medianFilterOrder=0;
tic
%[pitch, pitch2, pfMat]=ptByDpOverPf(y, fs, nbits, ptOpt, 1);
[pitch1, clarity1, pfMat1]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
fprintf('Time = %g seconds\n', toc);
time=(1:length(y))/fs;
frameTime=frame2sampleIndex(1:length(pitch1), ptOpt.frameSize, ptOpt.overlap)/fs;
%[pitch2, pitch3, pfMat]=ptByDpOverPf(y, fs, nbits, ptOpt);

subplot(4,1,1); plot(time, y); title(strPurify4label(waveFile));

subplot(4,1,2); plot(frameTime, pv, 'marker', 'o', 'color', 'g');
line(frameTime, pitch1, 'marker', '.', 'color', 'b');
%line(frameTime, pitch2, 'marker', '.', 'color', 'r');
%line(frameTime, pitch3, 'marker', '.', 'color', 'c');
title(sprintf('PT using ptByDpOverPfMex, with pfWeight=%d and indexDiffWeight=%d', ptOpt.pfWeight, ptOpt.indexDiffWeight));

ptOpt.useVolThreshold=1;
ptOpt.cutLeadingTrailingZero=0;
[pitch4, clarity2, pfMat2]=ptByPfMex(y, fs, nbits, ptOpt);
subplot(4,1,3); plot(frameTime, pv, 'marker', 'o', 'color', 'g');
line(frameTime, pitch4, 'marker', '.', 'color', 'b');
title('PT using ptByPfMex');

subplot(4,1,4)
plot(frameTime, clarity1, '.-', frameTime, clarity2, '.-');
line([min(frameTime), max(frameTime)], 0.5*max(clarity1)*[1 1], 'color', 'r');
set(gca, 'ylim', [0 1]);
 
return

pfType=1; ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.dpUseLocalOptim=1;
ptOpt.pfWeight=1;
weights=[0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
weights=0:3000:30000;
allPitch=[];
time=(1:length(y))/fs;

legendStr={};

subplot(2,1,1);
plot(time, y);  title(strPurify4label(sprintf('Waveform of %s', waveFile)));
set(gca, 'xlim', [-inf inf]);
subplot(2,1,2);
%plot(allPitch); legend(legendStr);
title(strPurify4label(sprintf('Pitch contours of %s', waveFile)));
box on;

tic
for i=1:length(weights)
	ptOpt.indexDiffWeight=weights(i);
	[pitch, clarity]=pitchTracking(waveFile, ptOpt);
	pitch(pitch==0)=nan;
	allPitch(:,i)=pitch(:)-2*(i-1);
	frameTime=frame2sampleIndex(1:length(pitch), ptOpt.frameSize, ptOpt.overlap)/fs;
	line(frameTime, allPitch(:,i), 'lineStyle', getLineStyle(i));
	index=find(~isnan(pitch));
	text(frameTime(index(end)), allPitch(index(end),i), ['  \theta=', int2str(weights(i))]); 
	
	legendStr{end+1}=sprintf('ptOpt.indexDiffWeight=%d\n', weights(i));
end
fprintf('Average time = %g seconds\n', toc/length(weights));
set(gca, 'xlim', [-inf inf]);
xlabel('Time (seconds)');
ylabel('Pitch (Semitones)');

tic

return

for i=1:2
	for j=1:2
		ptOpt.dpMethod=i;
		ptOpt.dpPosDiffOrder=j;
	%	tic, pitch2 = ptByDpOverPfMex(y, fs, nbits, ptOpt.frameSize, ptOpt.overlap); time2=toc;
	%	tic, pitch2 = ptByDpOverPfMex(y, fs, nbits, ptOpt.frameSize, ptOpt.overlap, ptOpt.smdfWeight, ptOpt.indexDiffWeight); time2=toc;
	%	tic, pitch2 = ptByDpOverPfMex(y, fs, nbits, ptOpt.frameSize, ptOpt.overlap, ptOpt.smdfWeight, ptOpt.indexDiffWeight, ptOpt.dpMethod, ptOpt.dpPosDiffOrder); time2=toc;
		tic, pitch2 = ptByDpOverPfMex(y, fs, nbits, ptOpt.frameSize, ptOpt.overlap, ptOpt.smdfWeight, ptOpt.indexDiffWeight, ptOpt.dpMethod, ptOpt.dpPosDiffOrder, ptOpt.dpUseLocalOptim); time2=toc;
		tic, [pitch1, volume1, volTh1, amdfMat1, A1] = wave2pitchByAmdfDp(y, fs, nbits, ptOpt); time1=toc;
		maxDiff=max(abs(pitch1-pitch2));
		fprintf('ptOpt.dpMethod=%d, ptOpt.dpPosDiffOrder=%d, m-time=%g, mex-time=%g, Max. diff in pitch=%g\n', i, j, time1, time2, maxDiff);
		
		if maxDiff>1
			subplot(2,1,1);
			plot(1:length(pitch1), pitch1, '.-', 1:length(pitch2), pitch2, '.-'); title('Identified pitch');
			legend('pitch by m', 'pitch by mex');
			subplot(2,1,2);
			plot(1:length(pitch1), abs(pitch1-pitch2), '.-'); title('Max. difference in pitch');
			return;
		end
	end
end