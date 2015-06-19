function optPath=wave2suv(waveObj, pdOpt, gmmSet, transProb, plotOpt, inputName)

if nargin<1, waveObj='\dataset\childSong4public\MIR-QBSH-corpus\waveFile\year2007/person00001/00011.wav'; end
if nargin<2, pdOpt=pdOptSet; end
if nargin<5, plotOpt=0; end
if nargin<6, inputName={'a', 'b'}; end

if ischar(waveObj),	waveObj=waveFile2obj(waveObj); end		% waveObj is actual the wave file name

[inData, outData, inputName, annotation]=wave2feature(waveObj, pdOpt);
%load bestInputIndex.mat
%inData=inData(bestInputIndex, :);
%inputName=inputName(bestInputIndex);
if pdOpt.useInputNormalize==1
	error('This is not implemneted yet!\n');
end

[feaDim, frameNum]=size(inData);
[optPath, stateProb, hmmTable]=hmm4suv(inData, gmmSet, transProb);

if plotOpt
	sampleTime=(1:length(waveObj.signal))'/waveObj.fs;
	frameTime=frame2sampleIndex(1:frameNum, pdOpt.frameSize, pdOpt.overlap)';
	subplot(5,1,1); plot(sampleTime, waveObj.signal);
	subplot(5,1,2); plot(frameTime, inData', '.-');  set(gca, 'xlim', [-inf inf]); legend(inputName);
	subplot(5,1,3); plot(frameTime, stateProb, '.-'); legend('SU', 'V'); set(gca, 'xlim', [-inf inf]);
	pitched=1+(stateProb(:,2)-stateProb(:,1)>0);
	subplot(5,1,4); plot(frameTime, outData, 'bo-', frameTime, optPath, 'k.-', frameTime, pitched, 'ro-'); legend('GT', 'Predcited', 'StateProb comp'); set(gca, 'xlim', [-inf inf]);
	%subplot(6,1,5); imagesc(stateProb'); axis xy
	subplot(5,1,5); imagesc(hmmTable'); axis xy
end
%keyboard
