% Collect the data for su/v detection

%	Roger Jang, 20060327, 20070417

pdOpt=pdOptSet;
fprintf('Collecting wave data and features from "%s"...\n', pdOpt.waveDir);
waveData=recursiveFileList(pdOpt.waveDir, 'wav');	% Collect all wave files
for i=1:length(waveData)
	fprintf('%d/%d ==> %s\n', i, length(waveData), waveData(i).name);
	waveFile=waveData(i).path;
	pvFile=strrep(waveFile, '.wav', '.pv');
	if ~exist(pvFile)
		fprintf('\tNo .pv file found for %s\n', waveFile);
		continue;
	end
	plotOpt=0;
	[waveData(i).inData, waveData(i).outData, inputName, waveData(i).annotation]=wave2feature(waveFile, pdOpt, plotOpt);
end

% Save the collect data to DS.mat
DS.inputName=inputName;
DS.outputName={'Unpitched', 'pitched'};
DS.input=[waveData.inData];
DS.output=[waveData.outData];
DS.annotation=[waveData.annotation];		% Pay attention: This is [], not {}!
fprintf('Save DS to DS.mat ...\n');
save DS.mat DS
