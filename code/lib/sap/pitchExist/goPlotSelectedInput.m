% 2D scatter plot

close all;
addMyPath;
load DS.mat
load bestInputIndex.mat
fprintf('bestInputIndex=%s\n', mat2str(bestInputIndex));
if length(bestInputIndex)~=2
	fprintf('No way to print scatter data of dim %d!\n', length(bestInputIndex));
	return;
end
DS.input=DS.input(bestInputIndex, :);			% Take the selected features
DS.inputName=DS.inputName(bestInputIndex);

figure; dcprDataPlot(DS, 'Original data');	% Original data plot
set(gcf, 'name', 'Original Data'); title('Original Data');
DS.input=inputNormalize(DS.input);
figure; dcprDataPlot(DS, 'Normalized data');	% Data plot after normalization
set(gcf, 'name', 'Normalized Data'); title('Normalized Data');