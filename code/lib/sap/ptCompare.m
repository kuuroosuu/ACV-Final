function [pitch1, pitch2]=ptCompare(waveFile)
%ptCompare: Comparison of pitch tracking methods

if nargin<1, waveFile='10LittleIndians.wav'; end

pvFile=  [waveFile(1:end-4), '.pv'];
pv=asciiRead(pvFile); 

[y, fs, nbits]=wavReadInt(waveFile);
wObj=waveFile2obj(waveFile);

pfType=1; ptOpt=ptOptSet(fs, nbits, pfType);
ptOpt.useVolThreshold=1;
ptOpt.cutLeadingTrailingZero=0;

time=(1:length(y))/fs;
subplot(4,1,1); plot(time, y); title(strPurify4label(waveFile));

[pitch1, clarity1, pfMat1]=ptByDpOverPfMex(y, fs, nbits, ptOpt);
frameTime=frame2sampleIndex(1:length(pitch1), ptOpt.frameSize, ptOpt.overlap)/fs;
subplot(4,1,2); plot(frameTime, pv, 'marker', 'o', 'color', 'g');
line(frameTime, pitch1, 'marker', '.', 'color', 'b');
title(sprintf('PT using ptByDpOverPfMex, with pfWeight=%d and indexDiffWeight=%d', ptOpt.pfWeight, ptOpt.indexDiffWeight));

%[pitch2, clarity2, pfMat2]=ptByPfMex(y, fs, nbits, ptOpt);
[pitch2, clarity2]=ptByPf(wObj, ptOpt);
subplot(4,1,3); plot(frameTime, pv, 'marker', 'o', 'color', 'g');
line(frameTime, pitch2, 'marker', '.', 'color', 'b');
title('PT using ptByPf');

subplot(4,1,4); plot(frameTime, clarity1, '.-', frameTime, clarity2, '.-');
title('Clarity using ptByDpOverPfMex & ptByPfMex');

sound(y/(2^nbits/2), fs);

% Create buttons
pitchObj1.frameRate=fs/(ptOpt.frameSize-ptOpt.overlap); pitchObj1.signal=pitch1;
pitchObj2=pitchObj1; pitchObj2.signal=pitch2;
buttonH=wavePitchPlayButton(wObj, pitchObj1, pitchObj2);