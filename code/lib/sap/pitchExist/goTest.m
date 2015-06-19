% Test the classifier for SU/V detection

addMyPath;			% Add required toolboxes to the search path
waveDir='testDataByByway';	% Directory for wave files for testing

waveData=recursiveFileList(waveDir, 'wav');			% Get all wave files under waveDir
load classifierParam.mat					% load bestInputIndex, qcParam, mu, sigma
[recogRate, DS, DS2]=getRecogRate(waveData, qcParam, bestInputIndex, mu, sigma);	% Perform feature extraction and recognition
fprintf('Recognition rate = %f%%\n', recogRate*100);

if length(bestInputIndex)==2
	decBoundaryPlot;
end