function pdOpt=pdOptSet
% pdOpt: Options for pitch detection (to determine if a frame has pitch or not)
%
%	Usage:
%		pdOpt=pdOptSet;
%
%	Description:
%		pdOpt=pdOptSet returns the default options for pitch detection.
%
%	Example:
%		pdOpt=pdOptSet

%	Category: Options for pitch detection
%	Roger Jang, 20130114

%% === Make sure you change the follow paths to where you store the following toolboxes
homeDir='d:/users/jang';
addpath([homeDir, '/matlab/toolbox/utility']);
addpath([homeDir, '/matlab/toolbox/sap']);
addpath([homeDir, '/matlab/toolbox/machineLearning']);
%% === Function for feature extraction
vdOpt.feaExtractFun=@pdFeaExtractFromFile;
%% === Folder for wave files
pdOpt.waveDir='\dataset\childSong4public\MIR-QBSH-corpus\waveFile\year2007\person00001';
%pdOpt.waveDir='\dataset\childSong4public\MIR-QBSH-corpus\waveFile\year2007';
%% === Parameters for PD
pdOpt.fs=8000;
pdOpt.frameSize=256;
pdOpt.overlap=0;
pdOpt.useInputNormalize=0;
pdOpt.frameZeroMeanOrder=3;