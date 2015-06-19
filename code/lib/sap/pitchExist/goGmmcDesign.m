% Classifier design based on GMM

addMyPath;
pdOpt=pdOptSet;

% ====== Load DS.mat
fprintf('Loading DS.mat...\n');
load DS.mat
load bestInputIndex.mat
%bestInputIndex=[1 2 3];		% This is for HMM
DS2=DS;								% DS2 is DS after input selection and normalization
DS2.input=DS2.input(bestInputIndex, :);				% Use the selected features
DS2.inputName=DS2.inputName(bestInputIndex);			% Update the input names based on the selected features
if pdOpt.useInputNormalize
	[DS2.input, mu, sigma]=inputNormalize(DS2.input);	% Input normalization
end
% ====== Down sample the training data
%fprintf('Down sampling DS2...\n');
%DS2.input=DS2.input(:, 1:10:end);
%DS2.output=DS2.output(:, 1:10:end);

% ====== GMM classifier
cvData=cvDataGen(DS2, 2, 'full');	% 2-fold cross-validation
TS=cvData(1).TS;
VS=cvData(1).VS;
gmmcOpt=gmmcOptSet;
gmmcOpt.config.gaussianNum=2:2:50;
%gmmcOpt.config.gaussianNum=[2, 4, 8, 16, 32, 64, 128];
gmmcOpt.config.covType=1;
showPlot=1;
figure;
[gmmData, recogRate1, recogRate2]=gmmcGaussianNumEstimate(TS, VS, gmmcOpt, showPlot);
[maxRr, index]=max(recogRate2);
fprintf('Saving gmmData...\n');
save gmmData gmmData index recogRate1 recogRate2
