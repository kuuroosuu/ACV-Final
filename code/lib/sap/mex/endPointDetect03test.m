% endPointDetect03.m 和 endPointDetect03Mex.dll 的同步測試

mex endpointDetect03Mex.cpp d:/users/jang/c/lib/audio.cpp d:/users/jang/c/lib/utility.cpp

addpath d:\users\jang\matlab\toolbox\audioProcessing
addpath d:\users\jang\matlab\toolbox\utility

% 單一檔案測試
%waveFile='D:\users\jang\application\pmpWaveFile\1.wav';
%[y, fs, nbits]=wavReadInt(waveFile);
%[endPoint, volume, zcr, soundSegment] = endPointDetect03(y, fs, nbits, 1);
%endPointIndex = endPointDetect03mex(y, fs);

% 多個檔案測試
waveFiles=recursiveFileList('D:\users\jang\application\pmpWaveFile', 'wav');
waveNum=length(waveFiles);
for i=1:waveNum
	waveFile=waveFiles(i).path;
	fprintf('%d/%d: %s\n', i, waveNum, waveFile);
	[y, fs, nbits]=wavReadInt(waveFile);
	[endPoint, volume, zcr, soundSegment]=endPointDetect03(y, fs, nbits);
	temp=[[soundSegment.beginFrame]; [soundSegment.endFrame]];
	epIndex1=temp(:)';
	epIndex2 = endPointDetect03mex(y, fs);
	same=isequal(epIndex1, epIndex2);
	if same~=1
		fprintf('Not equal!\n');
		break;
	end
end