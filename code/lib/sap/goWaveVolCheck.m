addpath /users/jang/matlab/toolbox/utility
addpath /users/jang/matlab/toolbox/dcpr


waveDir='D:\users\jang\temp\wave';
waveData=recursiveFileList(waveDir, 'wav');

for i=1:length(waveData)
	file=waveData(i).path;
	flag=waveVolumeCheck(file, [], 1);
	fprintf('%d/%d: file=%s, flag=%d\n', i, length(waveData), file, flag);
	fprintf('\tPress any key to continue...'); pause; fprintf('\n');
end