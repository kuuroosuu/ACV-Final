function [recogRate, DS, DS2]=suvPerfEval(pvData, suvOpt)
% suvPerfEval: Performance evaluation of SU/V detection

%	Roger Jang, 20070417

% Feature extraction from wave files
fprintf('Feature extraction...\n');
errorIndex=[];
for i=1:length(pvData)
	pvFile=pvData(i).path;
	waveFile = strrep(pvFile, '.pv', '.wav');
	fprintf('%d/%d ==> Collect features from %s\n', i, length(pvData), waveFile);
	[pvData(i).inData, pvData(i).outData, inputName, pvData(i).annotation]=feval(suvOpt.feaExtractFun, waveFile);
	if size(pvData(i).inData,2)~=size(pvData(i).outData,2)
		errorIndex=[errorIndex, i];
		fprintf('Error in reading %s (size(pvData(%d).inData,2)=%d, size(pvData(%d).outData,2)=%d)\n', pvData(i).path, i, size(pvData(i).inData,2), i, size(pvData(i).outData,2));
	end
end
fprintf('Removing %d error PV files...\n', length(errorIndex));
pvData(errorIndex)=[];

% Prepare DS
DS.inputName=inputName;
DS.input=[pvData.inData];
DS.output=double(cat(2, pvData.outData))+1;		% Cascade the desired output
DS.annotation=[pvData.annotation];			% Pay attention: This is [], not {}!

% Compute the recognition rate
fprintf('Compute recognition rate...\n');
[computedClass, DS2]=feval(suvOpt.mainClassifier, DS, suvOpt);			% Perform classification
recogRate=sum(computedClass==DS.output)/length(DS.output);	% Compute recognition rate