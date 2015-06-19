fprintf('Compiling wave2mfccMex.cpp...\n');
mex wave2mfccMex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/fft.cpp d:/users/jang/c/lib/myMfcc.cpp d:/users/jang/c/lib/utility.cpp

addpath d:/users/jang/matlab/toolbox/utility
addpath d:/users/jang/matlab/toolbox/audioProcessing

waveDir='C:\temp\【鈦映科技】語者確認錄音';
waveData=recursiveFileList(waveDir, 'wav');
caseNum = 5;

for i=1:caseNum
	waveFile = waveData(i).path;
	fprintf('%d/%g: %s\n', i, caseNum, waveFile);
	[y, fs]=wavread(waveFile);
	mp=mfccPrmSet(fs);

	% ====== wav2ftr1 is a M file that is original feature_extrac.m
	tic; mfcc1 = wave2mfcc(y, fs, mp); time1=toc;
	load mfcc.mat
	fprintf('isequal(mfcc1, mfcc)=%d\n', isequal(mfcc1, mfcc));
	% ====== A single dll
	tic; mfcc2 = wave2mfccmex(y, fs, mp); time2=toc;
	
	fprintf('time1 (m) = %g, time2 (dll) = %g\n', time1, time2);
	
	if ~isequal(size(mfcc1), size(mfcc2))
		fprintf('Dimensions mismatch:\n');
		fprintf('\tsize(mfcc1)=%s\n', mat2str(size(mfcc1)));
		fprintf('\tsize(mfcc2)=%s\n', mat2str(size(mfcc2)));
	else
		fprintf('size(mfcc1)=size(mfcc2)=%s\n', mat2str(size(mfcc1)));
		fprintf('Diff. between mfcc1 and mfcc2 = %g\n', sum(sum(abs(mfcc1-mfcc2))));
	end
	fprintf('\n');
end