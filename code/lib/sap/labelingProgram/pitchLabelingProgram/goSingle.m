% Label the pitch of a single wave file

addpath d:/users/jang/matlab/toolbox/utility -end
addpath d:/users/jang/matlab/toolbox/sap -end
close all; clear all

waveFile='\users\jang\books\audioSignalProcessing\programmingContest\pitchTracking\exampleProgram\waveFile\rogerJang\old macdonald had a farm_unknown_0.wav';
au=myAudioRead(waveFile);
ptOpt=ptOptSet(au.fs, au.nbits);
ptOpt.waveFile=au.file;
ptOpt.targetPitchFile=[ptOpt.waveFile(1:end-3), 'pv'];
ptOpt.frameSize=round(ptOpt.frameDuration*au.fs/1000);
ptOpt.overlap=round(ptOpt.overlapDuration*au.fs/1000);
ptOpt.frame2pitchFcn='frame2pitch4labeling';
pitchLabel(au, ptOpt);