% epdByVol2.m 和 epdByVol2mex.dll 的同步測試
mex epdByVol2mex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

addpath d:\users\jang\matlab\toolbox\audioProcessing
addpath d:\users\jang\matlab\toolbox\utility

% 單一檔案測試
%waveFile='D:\users\jang\application\pmpWaveFile\1.wav';
%[y, fs, nbits]=wavReadInt(waveFile);
%[endPoint, volume, zcr, soundSegment] = epdByVol2(y, fs, nbits, 1);
%endPointIndex = epdByVol2mex(y, fs);

% 多個檔案測試
waveFiles=recursiveFileList('C:\dataSet\資策會-派斌\pmpWaveFile', 'wav');
waveNum=length(waveFiles);
for i=1:waveNum
	waveFile=waveFiles(i).path;
	fprintf('%d/%d: %s\n', i, waveNum, waveFile);
	[y, fs, nbits]=wavReadInt(waveFile);
	[endPoint, volume, zcr, soundSegment]=epdByVol2(y, fs, nbits, 1);
	temp=[[soundSegment.beginFrame]; [soundSegment.endFrame]];
	epIndex1=temp(:)';
	epIndex2 = epdByVol2mex(y, fs);
	same=isequal(epIndex1, epIndex2);
	if same~=1
		fprintf('Not equal!\n');
	else
		fprintf('Equal!\n');
	end
	fprintf('Press any key to continue...'); pause; fprintf('\n');
end