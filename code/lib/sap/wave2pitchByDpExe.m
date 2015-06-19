function [pitch, pitch2, volume]=wave2pitchByDpExe(waveFile, ptOpt, plotOpt)
% wave2pitch: Pitching tracking using DP over amdf matrix
%	Usage: [pitch, pitch2, volume]=wave2pitch(waveFile, ptOpt, plotOpt)
%		pitch: pitch vector of the whole wave file
%		pitch2: pitch vector after volume thresholding
%		ptOpt.indexDiffWeight and ptOpt.pfWeight are parameters for DP
%		This function required the use of wave2pitchVolume.exe.
%
%	For example:
%		waveFile='二十四橋明月夜.wav';
%		[y, fs, nbits]=wavread(waveFile);
%		ptOpt=ptOptSet(fs, nbits);
%		plotOpt=1;
%		[pitch, pitch2, volume]=wave2pitchByDpExe(waveFile, ptOpt, plotOpt);

%	Roger Jang, 20060218

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end
debug=0;

% ====== Resample the wave file if necessary
tempWaveFile=[tempname, '.wav'];
wave2wave(waveFile, tempWaveFile, 16000, 16);	% Convert to 16KHz, 16Bits

[y, fs, nbits]=wavReadInt(tempWaveFile);
% ====== 抓出整段音高及音量
[parentDir, junk, junk]=fileparts(which(mfilename));
exeDir=[parentDir, '\exe'];
executable=[exeDir, '\wave2pitchVolume.exe'];
if ~exist(executable)
	msg=sprintf('Cannot find %s!\n', executable);
	error(msg);
end
pitchFile=[tempname, '.pv'];
pitchFile2=[tempname, '.pv2'];
volumeFile=[tempname, '.vol'];
useEpd=1;
if isempty(ptOpt.indexDiffWeight)
	dosCmd=sprintf('%s "%s" -1 -1 %d %d %s %s', executable, tempWaveFile, ptOpt.frameSize, ptOpt.overlap, pitchFile, volumeFile);
else
	dosCmd=sprintf('%s "%s" -1 -1 %d %d %s %s %g %g %d %s', executable, tempWaveFile, ptOpt.frameSize, ptOpt.overlap, pitchFile, volumeFile, ptOpt.indexDiffWeight, ptOpt.pfWeight, useEpd, pitchFile2);
end
%dosCmd
currDir=pwd; cd(exeDir);
[exeStatus, exeResult]=dos(dosCmd);
cd(currDir);
if debug
	fprintf('dosCmd=%s\n', dosCmd);
	fprintf('exeResult=%s\n', exeResult);
end

if exist(pitchFile)~=2 | exist(pitchFile2)~=2 | exist(volumeFile)~=2
	fprintf('Some output files are missing!\n');
	fprintf('Please check the result of running the following command within DOS:\n');
	fprintf('%s\n', dosCmd);
	return;
end

pitch=asciiRead(pitchFile); delete(pitchFile);
pitch2=asciiRead(pitchFile2); delete(pitchFile2);
volume=asciiRead(volumeFile); delete(volumeFile);
frameNum=length(pitch);
frameTime=frame2sampleIndex(1:frameNum, ptOpt.frameSize, ptOpt.overlap)/fs;
% ====== Plotting
if plotOpt
	plotNum=4;
	subplot(plotNum, 1, 1);
	time=(1:length(y))/fs;
	plot((1:length(y))/fs, y);
	axis([min(time) max(time) -2^nbits/2 2^nbits/2]); title(strPurify(waveFile)); grid on

	subplot(plotNum, 1, 2);
	plot(frameTime, volume, '.-');
	set(gca, 'xlim', [min(frameTime) max(frameTime)]); title('Volume'); grid on

	subplot(plotNum,1, 3);
	pfMat = wave2pfMat(y, ptOpt.frameSize, ptOpt.overlap, ptOpt.pfLen, 0, 3);	% 0 for AMDF, 1 for ACF
	pp=round(fs./pitch2freq(pitch));
	pfPlot(pfMat, frameTime, pp); ylabel('PF matrix');

	subplot(plotNum, 1, 4);
	pitchA=pitch; pitchA(pitchA==0)=nan;
	pitchB=pitch2; pitchB(pitchB==0)=nan;
	plot(frameTime, pitchA, '.-'); axis([min(frameTime) max(frameTime) 35 100]); title('Pitch'); grid on
	line(frameTime, pitchB, 'color', 'r', 'marker', 'o', 'linestyle', 'none');
%	pitch(pitch==0)=nan;
%	plot(frameTime, pitch, '.-'); axis([min(frameTime) max(frameTime) 35 100]); title('Pitch'); grid on
	% ====== Buttons for playback
	waveObj=waveFile2obj(waveFile);
%	waveObj.signal=y-mean(y); waveObj.fs=fs; waveObj.nbits=nbits;
	pitchObj.signal=pitchA; pitchObj.frameRate=fs/(ptOpt.frameSize-ptOpt.overlap);
	pitchObj2.signal=pitchB; pitchObj2.frameRate=fs/(ptOpt.frameSize-ptOpt.overlap);
	buttonH=wavePitchPlayButton(waveObj, pitchObj, pitchObj2);
	legend('Pitch: whole', 'Pitch2: segmented');
end

% ====== Selfdemo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
waveFile='dan_shi_lung_cheng_fei_qiang_zai.wav';
%waveFile='file56.wav';
%waveFile='file19.wav';
[y, fs, nbits]=wavread(waveFile);
ptOpt=ptOptSet(fs, nbits);	% Why pfType=0 gives wrong result? 
plotOpt=1;
[pitch, pitch2, volume]=wave2pitchByDpExe(waveFile, ptOpt, plotOpt);

return

% Plot the influence of ptOpt.indexDiffWeight on pitch tracking
weights=[1, 10, 100, 1000];
allPitch=[];
legendStr={};
for i=1:length(weights)
	ptOpt.indexDiffWeight=weights(i);
	thePitch=wave2pitchByDpExe(waveFile, ptOpt, 0)+i/2;
	allPitch(:,i)=thePitch(:);
	legendStr{end+1}=sprintf('ptOpt.indexDiffWeight=%d\n', weights(i));
end
figure; plot(allPitch);
legend(legendStr);