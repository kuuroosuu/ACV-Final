% Design a quadratic classifier
pdOpt=pdOptSet;
%% ====== Load DS.mat
fprintf('Loading DS.mat...\n');
load DS.mat
load bestInputIndex.mat
DS2=DS;								% DS2 is DS after input selection and normalization
DS2.input=DS2.input(bestInputIndex, :);				% Use the selected features
DS2.inputName=DS2.inputName(bestInputIndex);			% Update the input names based on the selected features
if pdOpt.useInputNormalize
	[DS2.input, mu, sigma]=inputNormalize(DS2.input);	% Input normalization
end
priors=dsClassSize(DS2);					% Take data count of each class as the class probability
%% ====== Down sample the training data
%fprintf('Down sampling DS2...\n');
%DS2.input=DS2.input(:, 1:10:end);
%DS2.output=DS2.output(:, 1:10:end);
%% ====== Qudratic classifier
[qcParam, logLike, recogRate]=qcTrain(DS2, [], 1);		% Classifier design based on quadratic classifiers. Also plot the data distribution
fprintf('Recognition rate = %f%%\n', recogRate*100);

fprintf('Saving classifier-related info to classifierParam.mat...\n');
if pdOpt.useInputNormalize
	save qcParam.mat qcParam bestInputIndex mu sigma
else
	save qcParam.mat qcParam bestInputIndex
end

if length(bestInputIndex)==2
	decBoundaryPlot;
end
