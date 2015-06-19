% 對每個 wave 檔案進行端點偵測，顯示結果，並將結果寫入特定目錄

waveDir='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機';
waveDir=strrep(waveDir, '\', '/');
if isempty(which('findcellstr')); addpath('d:/users/jang/matlab/toolbox/utility'); end
plotOpt=0;
waveFiles=dir([waveDir, '/*.wav']);
for i=1:length(waveFiles)
	fileName=[waveDir, '/', waveFiles(i).name];
	answer=file2answer(fileName);
	fprintf('\n%g/%g: Hit return to play %s (%s):', i, length(waveFiles), fileName, answer);
	if plotOpt, pause; end
	[y, fs, nbits]=wavread(fileName);
	output=endPoint2(y, fs, nbits, [], plotOpt);
	startIndex=output(1);
	endIndex=output(2);
	if plotOpt, subplot(3,1,1); title([fileName, ' (', answer, ')']); end
	newY=y(startIndex:endIndex);
	targetFile=[waveDir, '/cleaned/', waveFiles(i).name];
	wavwrite(newY, fs, nbits, targetFile);
end

%return

% 測試辨識率
currDir=pwd;
cd D:\users\jang\matlab\toolbox\asr\phoneRecog
cmd='scoring D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\cleaned';
dos(cmd);
cd(currDir);