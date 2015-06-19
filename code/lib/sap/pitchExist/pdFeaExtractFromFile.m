function [feature, output, inputName, annotation, DS]=pdFeaExtractFromFile(wObj, pdOpt, showPlot)
% pdFeaExtractFromFile: Wave to feature conversion for s/u/v detection
%
%	Usage:
%		[feature, output, inputName, annotation]=pdFeaExtractFromFile(wObj, showPlot)
%
%	Example:
%		waveFile='d:/dataSet/childSong/waveFile/2007-音訊處理與辨識/19461108任佳王民/十個印第安人_不詳_0.wav';
%		pdOpt=pdOptSet;
%		[feature, output, inputName, annotation]=pdFeaExtractFromFile(waveFile, pdOpt, 1);

%	Category: Feature extraction for pitch detection
%	Roger Jang, 20040910, 20070417

if nargin<1, selfdemo; return; end
if nargin<2, pdOpt=pdOptSet; end
if nargin<3, showPlot=0; end

if ischar(wObj), wObj=waveFile2obj(wObj); end		% wObj is actual the wave file name
wObj=waveFormatConvert(wObj, 8000, 8, 1);				% Format conversion to mono, 8KHz, 8bits

% ====== Frame blocking
frameSize=pdOpt.frameSize;
overlap=pdOpt.overlap;
frameMat=buffer2(wObj.signal, frameSize, overlap);
frameNum=size(frameMat, 2);
% ====== zero-mean for each frame
frameMat=frameZeroMean(frameMat, pdOpt.frameZeroMeanOrder);

feature=[];
inputName={};

% ====== Feature 1: volume
vol = frame2volume(frameMat);
volRatio=vecRatio(vol);
feature=[feature; volRatio];
inputName={inputName{:}, 'vol'};
% ====== Feature 1: zcr
%zcr = frame2zcr(frameMat, 2);
%zcrRatio=vecRatio(zcr);
%feature=[feature; zcrRatio];
%inputName={inputName{:}, 'zcr'};
% ====== Feature 2: clarity
clarity=frame2clarity(frameMat, wObj.fs, 'acf');
clarityRatio=vecRatio(clarity);
feature=[feature; clarityRatio];
inputName={inputName{:}, 'clarity'};
% ====== Feature 3: hod
ashod = frame2ashod(frameMat, 4);
ashod=vecRatio(ashod);
feature=[feature; ashod];
inputName={inputName{:}, 'hod'};

% ====== add annotation
for i=1:size(feature,2)
	annoStr=sprintf('%s\n%s', wObj.file, int2str(i));
	annoStr=strrep(annoStr, '\', '/');
	annoStr=strrep(annoStr, '_', '\_');
	annotation{i}=annoStr;
end

DS.inputName=inputName;
DS.input=feature;
DS.annotation=annotation;

% ====== Read human-labeled pitch file
output=[];
pvFile=strrep(wObj.file, '.wav', '.pv');
if exist(pvFile)==2
	targetPitch=asciiRead(pvFile);
	if length(targetPitch)>frameNum, targetPitch=targetPitch(1:frameNum); end	% Due to the difference between buffer.m (used before) and buffer2.m (used now)
	output=targetPitch>0;
	output=output+1;	% {0,1} ===> {1,2}
	output=output(:)';
end
DS.output=output;

if showPlot
	plotTitle=strrep(wObj.file, '\', '/'); plotTitle=strrep(plotTitle, '_', '\_');
	DS.input=feature(1:2, :);
	DS.output=output;
	DS.inputName=inputName;
	DS.annotation=annotation;
	subplot(211); dsScatterPlot(DS); title(['Raw data: ', plotTitle]);
	feature2=inputNormalize(feature(1:2, :));
	DS.input=feature2(1:2, :);
	DS.output=output;
	subplot(212); dsScatterPlot(DS); title(['Normalized data: ', plotTitle]);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);