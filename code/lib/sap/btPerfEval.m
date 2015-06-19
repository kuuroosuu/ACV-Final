function [meanFMeasure, waveData, totalInvalidCount]=btPerfEval(waveData, btOpt)
% Roger Jang, 20110811

if nargin<2, btOpt=btOptSet; end

waveNum=length(waveData);
% ====== Beat tracking
%fprintf('\tBeat tracking of %d wave files from "%s"...\n', waveNum, btOpt.waveDir);
for i=1:waveNum
%	fprintf('\t%d/%d ==> waveFile=%s, ', i, waveNum, waveData(i).file);
	tic;
	if ~isfield(waveData, 'osc'), waveData(i).osc=wave2osc(waveData(i).wObj, btOpt.oscOpt); end
	switch(btOpt.type)
		case 'time-varying'
			if ~isfield(waveData, 'tempoCurve'), waveData(i).tempoCurve=tempogram2tempoCurve(waveData(i).tempogram, btOpt); end
			if ~isfield(waveData, 'tempogram'), waveData(i).tempogram=osc2tempogram(waveData(i).osc, btOpt); end
			waveData(i).cBeat=tempoCurve2beat(waveData(i).tempoCurve, waveData(i).osc, waveData(i).tempogram, btOpt);
		case 'constant'
			waveData(i).cBeat=feval(btOpt.btFcn, waveData(i), btOpt);
	end
	% === Performance evaluation
	for j=1:length(waveData(i).gtBeat)
		waveData(i).fMeasure(j)=simSequence(waveData(i).cBeat, waveData(i).gtBeat{j}, btOpt.tolerance);
	end
	waveData(i).meanFMeasure=nanmean(waveData(i).fMeasure);
	waveData(i).invalidCount=sum(isnan(waveData(i).fMeasure));
%	fprintf('mean F-measure=%.2f, time=%.2f sec\n', waveData(i).meanFMeasure, toc);
end
meanFMeasure=mean([waveData.meanFMeasure]);
totalInvalidCount=sum([waveData.invalidCount]);
%fprintf('Mean F-measure=%.2f\n', meanFMeasure);

% ====== Save the collect data
%fprintf('Save waveData to %s...\n', btOpt.waveDataFile);
%save(btOpt.waveDataFile, 'waveData');
