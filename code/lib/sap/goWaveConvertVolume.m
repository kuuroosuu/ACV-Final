% Convert wave files to 16k, 16b, mono
% Note the original file will be OVERWRITTEN!!!

addpath /users/jang/matlab/toolbox/utility
addpath /users/jang/matlab/toolbox/audioProcessing

waveDir='C:\temp\Sophia';
fprintf('Collect wave file info from %s...\n', waveDir);
waveData=recursiveFileList(waveDir, 'wav');
waveNum=length(waveData);

for i=1:waveNum
	source=waveData(i).path;
	fprintf('%d/%d: %s\n', i, waveNum, source);
	[y, fs, nbits]=wavread(source);
	if fs==16000 & nbits==16 & size(y,2)==1
		continue;
	end
	yMax=max(abs(y));
	y=y/yMax;
	target=source;		% Overwrite to the original file
	try
		wavwrite(y, fs, nbits, target);
	catch
		fprintf('%s: Error!\n', source);
	end
end