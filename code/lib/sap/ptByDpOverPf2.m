function [pitch, pitch2, pfMat]=ptByDpOverPf(y, fs, nbits, ptOpt, epdPrm, plotOpt)
% ptByDpOverPf: Wave to pitch conversion using dynamic programming
%	Usage: [pitch, pitch2, pfMat]=ptByDpOverPf(y, fs, nbits, ptOpt, epdPrm, plotOpt)
%		pitch: segmented pitch (if ptOpt.useEpd is 1)
%		pitch2: whole pitch
%		pfMat: pitch function matrix for DP
%
%	For example:
%		% ====== For speech input
%		waveFile='what_movies_have_you_seen_recently.wav';
%		[y, fs, nbits]=wavReadInt(waveFile);
%		ptOpt=ptOptSet(fs, nbits);
%		ptOpt.useEpd=1;
%		plotOpt=1;
%		pitch=ptByDpOverPf(y, fs, nbits, ptOpt, [], plotOpt);
%		% ====== For singing input
%		waveFile='yankee_doodle.wav';
%		[y, fs, nbits]=wavReadInt(waveFile);
%		ptOpt=ptOptSet(fs, nbits);
%		plotOpt=1;
%		pitch=ptByDpOverPf(y, fs, nbits, ptOpt, [], plotOpt);
%
%	Roger Jang, 20070824

if nargin<1, selfdemo; return; end
if nargin<2, fs=8000; end
if nargin<3, nbits=8; end
if nargin<4 | isempty(ptOpt), ptOpt=ptOptSet(fs, nbits); end
if nargin<5 | isempty(epdPrm), epdPrm=epdPrmSet(fs); end
if nargin<6, plotOpt=0; end

%{
fs2=8000;
if fs~=fs2
	fprintf('Convert fs from %d to %d...\n', fs, fs2);
	y=resample(y, fs2, fs);
end
fs=fs2;

nbits2=8;
if nbits~=nbits2
%	fprintf('Convert nbits from %d to %d...\n', nbits, nbits2);
	y=y/256;
end
nbits=nbits2;
%}

if ptOpt.useEpd
	if ptOpt.useEpdBeforePt	% 先做 EPD，再分段進行 PT
		% EPD
		[out1, out2, segment] = epdByVolHod(y, fs, nbits, epdPrm);
		frameNum1=floor((length(y)-epdPrm.overlap)/(epdPrm.frameSize-epdPrm.overlap));
		segmentNum=length(segment);
		% PT
		frameNum2=floor((length(y)-ptOpt.overlap)/(ptOpt.frameSize-ptOpt.overlap));
		pitch=zeros(1, frameNum2);
		clarity=zeros(1, frameNum2);
		pfMat=nan*ones(ptOpt.pfLen, frameNum2);
		for i=1:segmentNum
			index1=frame2sampleIndex(segment(i).beginFrame, epdPrm.frameSize, epdPrm.overlap);
			index2=frame2sampleIndex(segment(i).endFrame,   epdPrm.frameSize, epdPrm.overlap);
			[thePitch, theClarity, thePfMat]=ptByDpOverPfMex(y(index1:index2), fs, nbits, ptOpt);
			frameIndex=sample2frameIndex(index1, ptOpt.frameSize, ptOpt.overlap);
			pitch(frameIndex:frameIndex+length(thePitch)-1)=thePitch;
			clarity(frameIndex:frameIndex+length(thePitch)-1)=theClarity;
			pfMat(:, frameIndex:frameIndex+length(thePitch)-1)=thePfMat;
		end
		pitch2=pitch;		% In this case, pitch2 = pitch
	else % 先做 PT，再做 EPD 來判斷靜音
		% PT
		[pitch2, clarity2, pfMat]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
		% EPD
		[out1, out2, segment] = epdByVolHod(y, fs, nbits, epdPrm);
		frameNum1=floor((length(y)-epdPrm.overlap)/(epdPrm.frameSize-epdPrm.overlap));
		segmentNum=length(segment);
		pitch=0*pitch2;
		for i=1:segmentNum
			sId1=frame2sampleIndex(segment(i).beginFrame, epdPrm.frameSize, epdPrm.overlap);
			sId2=frame2sampleIndex(segment(i).endFrame,   epdPrm.frameSize, epdPrm.overlap);
			fId1=sample2frameIndex(sId1, ptOpt.frameSize, ptOpt.overlap);
			fId2=sample2frameIndex(sId2, ptOpt.frameSize, ptOpt.overlap);
			fId2=min(fId2, length(pitch));
			pitch(fId1:fId2)=pitch2(fId1:fId2);
		end
	end
else
	[pitch, clarity, pfMat]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
	pitch2=pitch;		% In this case, pitch2 = pitch
end
if ptOpt.useSmooth, pitch=pitchSmooth(pitch); end

% 找出 pitch 的 median，重新設定 ptOpt.minPitch, ptOpt.maxPitch, 再做一遍
if ptOpt.useRepeat
	pitch(pitch==0)=nan;
	medianValue=nanmedian(pitch);
	ptOpt.minPitch=medianValue-ptOpt.pitchRangeMax/2;
	ptOpt.maxPitch=medianValue+ptOpt.pitchRangeMax/2;
	if ptOpt.useEpd
		if ptOpt.useEpdBeforePt	% 先做 EPD，再分斷進行 PT
			% PT
			frameNum2=floor((length(y)-ptOpt.overlap)/(ptOpt.frameSize-ptOpt.overlap));
			pitch=zeros(1, frameNum2);
			pfMat=nan*ones(ptOpt.pfLen, frameNum2);
			for i=1:segmentNum
				index1=frame2sampleIndex(segment(i).beginFrame, epdPrm.frameSize, epdPrm.overlap);
				index2=frame2sampleIndex(segment(i).endFrame,   epdPrm.frameSize, epdPrm.overlap);
				[thePitch, theClarity, thePfMat]=ptByDpOverPfMex(y(index1:index2), fs, nbits, ptOpt);
				frameIndex=sample2frameIndex(index1, ptOpt.frameSize, ptOpt.overlap);
				pitch(frameIndex:frameIndex+length(thePitch)-1)=thePitch;
				clarity(frameIndex:frameIndex+length(thePitch)-1)=theClarity;
				pfMat(:, frameIndex:frameIndex+length(thePitch)-1)=thePfMat;
			end
			pitch2=pitch;
		else % 先做 PT，再做 EPD 來判斷靜音
			% PT
			[pitch2, clarity2, pfMat]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
			pitch=0*pitch2;
			for i=1:segmentNum
				sId1=frame2sampleIndex(segment(i).beginFrame, epdPrm.frameSize, epdPrm.overlap);
				sId2=frame2sampleIndex(segment(i).endFrame,   epdPrm.frameSize, epdPrm.overlap);
				fId1=sample2frameIndex(sId1, ptOpt.frameSize, ptOpt.overlap);
				fId2=sample2frameIndex(sId2, ptOpt.frameSize, ptOpt.overlap);
				fId2=min(fId2, length(pitch));
				pitch(fId1:fId2)=pitch2(fId1:fId2);
			end
		end
	else
		[pitch, clarity, pfMat]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
		pitch2=pitch;
	end
	if ptOpt.useSmooth, pitch=pitchSmooth(pitch); end
end

if plotOpt
	subplot(3,1,1);
	time=(1:length(y))/fs;
	plot(time, y);
	
	if ptOpt.useEpd
		frameTime1=frame2sampleIndex(1:frameNum1, epdPrm.frameSize, epdPrm.overlap)/fs;
		for i=1:length(segment)
			line(frameTime1(segment(i).beginFrame)*[1 1], 2^nbits/2*[-1, 1], 'color', 'm');
			line(frameTime1(segment(i).endFrame  )*[1 1], 2^nbits/2*[-1, 1], 'color', 'g');
		end
		ylabel('Amplitude'); title('Waveform');
	end
	axis([min(time) max(time) -2^nbits/2, 2^nbits/2]);
	frameNum2=floor((length(y)-ptOpt.overlap)/(ptOpt.frameSize-ptOpt.overlap));
	frameTime2=frame2sampleIndex(1:frameNum2, ptOpt.frameSize, ptOpt.overlap)/fs;
	subplot(3,1,2);
	imagesc(pfMat); axis xy
	subplot(3,1,3);
	temp=pitch; temp(temp==0)=nan;
	temp2=pitch2; temp2(temp2==0)=nan;
	plot(frameTime2, temp2, 'go-', frameTime2, temp, 'k.');
	if ptOpt.useRepeat
		line([min(frameTime2), max(frameTime2)], medianValue*[1 1], 'color', 'm');
		line([min(frameTime2), max(frameTime2)], ptOpt.minPitch*[1 1], 'color', 'm');
		line([min(frameTime2), max(frameTime2)], ptOpt.maxPitch*[1 1], 'color', 'm');
	end
	axis([min(frameTime2) max(frameTime2) -inf inf]);
	xlabel('Time in seconds'); ylabel('Semitones'); title('Pitch'); grid on
	% ====== Buttons for playback
	waveObj.signal=y-mean(y); waveObj.fs=fs; waveObj.nbits=nbits;
	pitchObj.signal=pitch; pitchObj.frameRate=fs/(ptOpt.frameSize-ptOpt.overlap);
	buttonH=wavePitchPlayButton(waveObj, pitchObj);
end

% ====== Selfdemo
function selfdemo
waveFile='lately2.wav';		% 16 kHz, 16 bits
waveFile='但使龍城飛將在.wav';
waveFile='what_movies_have_you_seen_recently.wav';
%waveFile='yankee_doodle.wav';
%waveFile='twinkle_twinkle_little_star.wav';
[y, fs, nbits]=wavReadInt(waveFile);
pfType=1;	% 0 for AMDF, 1 for ACF
ptOpt=ptOptSet(fs, nbits, pfType);
epdPrm=epdPrmSet(fs);
plotOpt=1;
[pitch, pitch2, pfMat]=ptByDpOverPf2(y, fs, nbits, ptOpt, epdPrm, plotOpt);