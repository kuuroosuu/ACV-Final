
waveFile='\dataset\childSong4public\MIR-QBSH-corpus\waveFile\year2007/person00001/00011.wav';
[y, fs, nbits]=wavReadInt(waveFile);

pdOpt=pdOptSet;
%load bestInputIndex.mat
[inData, outData, inputName, annotation]=wave2feature(waveFile, pdOpt);
%inData=inData(bestInputIndex, :);
%inputName=inputName(bestInputIndex);
if pdOpt.useInputNormalize==1
	error('This is not implemneted yet!\n');
end

[feaDim, frameNum]=size(inData);

load gmmData.mat
load transProb.mat
transLogProb=log(transProb);
%transLogProb=zeros(size(transProb));
gmmcModel=gmmData(index);
%[optPath, stateProb]=hmm4suv(inData, gmmcModel, transLogProb);
suv=wave2suv(waveFile, pdOpt, gmmcModel, transProb, 1, inputName);

% ====== GMM recog. rate
computedClassIndex = gmmClassifierEval(inData, gmmcModel);
rr=sum(outData==computedClassIndex)/length(outData);
fprintf('rr=%.2f%%\n', rr*100);
