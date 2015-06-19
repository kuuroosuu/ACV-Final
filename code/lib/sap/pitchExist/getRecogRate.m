function [recogRate, DS, DS2]=getRecogRate(waveData, qcParam, bestInputIndex, mu, sigma)
% getRecogRate: Compute the recognition rate of SU/V detection

%	Roger Jang, 20070417

% Feature extraction from wave files
fprintf('Feature extraction...\n');
for i=1:length(waveData)
	waveFile=waveData(i).path;
	plotOpt=0;
	[waveData(i).inData, waveData(i).outData, inputName, waveData(i).annotation]=wave2feature(waveFile, plotOpt);
end

% Prepare DS
DS.inputName=inputName;
DS.input=[waveData.inData];
DS.output=double(cat(2, waveData.outData));		% Cascade the desired output
DS.annotation=[waveData.annotation];			% Pay attention: This is [], not {}!
dataNum=size(DS.input, 2);

% Prepare DS2 (DS after input selection and normalization)
DS2=DS;
DS2.input=DS2.input(bestInputIndex, :);				% Use the selected features
DS2.inputName=DS2.inputName(bestInputIndex);			% Update the input names based on the selected features
DS2.input=inputNormalize(DS2.input, 0, mu, sigma);			% Input normalization based on the given mu and sigma

% Compute the recognition rate
fprintf('Compute recognition rate...\n');
for i=1:dataNum
%	fprintf('%d/%d\n', i, dataNum);
	DS3.input=DS2.input(:,i);					% Get an entry
	computedClass(i)=myClassifier(DS3, qcParam);			% Perform classification
end
recogRate=sum(computedClass==DS.output)/length(DS.output);	% Compute recognition rate