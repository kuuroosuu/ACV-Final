vdOpt=vdOptSet;
% === Collect wave data
waveData=waveDataFeaCollect(vdOpt);
% === HMM train
hmmModel=vdHmmTrain(waveData, vdOpt, 1);
% === HMM test
waveFile='D:\dataset\vibrato\female\combined-female.wav';
cOutput=vdHmmEval(waveFile, vdOpt, hmmModel, 1);