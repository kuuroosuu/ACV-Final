% Input selection based on KNNR and LOO

addMyPath;
fprintf('Loading DS.mat...\n');
load DS.mat				% Load the dataset
% ====== Down sample the training data
%fprintf('Down sample DS.mat...\n');
%DS.input=DS.input(:, 1:100:end);
%DS.output=DS.output(:, 1:100:end);

%DS.input=inputNormalize(DS.input);	% Input (feature) normalization
inputNum=3;				% Select 2 features at most
plotOpt=1;				% Plot the result of input selection
classifier='knncLoo';			% Use knnrLoo.m for input selection
param.k=1;				% Parameters for the classifier
fprintf('inputName=%s\n', cell2str(DS.inputName));	% Print the input names
selectionMethod='exhaustive';		% Method for input selection
switch selectionMethod
	case 'sequential'
		% ====== Sequential forward selection
		[bestInputIndex, bestRecogRate, allSelectedInput, allRecogRate, elapsedTime] = inputSelectSequential(DS, inputNum, classifier, param, plotOpt);
	case 'exhaustive'
		% ====== Exhaustive search
		[bestInputIndex, bestRecogRate, allSelectedInput, allRecogRate, elapsedTime] = inputSelectExhaustive(DS, inputNum, classifier, param, plotOpt);
	otherwise
		disp('Unknown method!');
end

fprintf('Recognition rate = %f\n', bestRecogRate*100);
fprintf('bestInputIndex = %s\n', mat2str(bestInputIndex));
fprintf('bestInputIndex = %s\n', cell2str(DS.inputName(bestInputIndex)));

[allRecogRate, index]=sort(allRecogRate);
allSelectedInput=allSelectedInput(index);
figure; inputSelectPlot(allRecogRate*100, allSelectedInput, DS.inputName, mfilename);
fprintf('Save the indices of the selected features to bestInputIndex.mat...\n');
save bestInputIndex bestInputIndex