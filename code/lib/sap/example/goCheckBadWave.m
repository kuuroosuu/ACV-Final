% Check wrongly recognized wave files (listed in badWave.txt) in a given directory. 

directory='D:\users\jang\matlab\toolbox\asr\application\anYuan\waveData\簡國慶2';
directory='D:\users\jang\matlab\toolbox\asr\application\RMC\waveData\ivr辨識音檔-手機\cleaned';
if isempty(which('findcellstr')); addpath('d:/users/jang/matlab/toolbox/utility'); end
badWaveFile=[directory, '/badWave.txt'];
contents=fileRead(badWaveFile);
%directory=contents{1};
contents(1)=[];
for i=1:length(contents)
	items=split(contents{i}, ' ');
	waveFile=[directory, '/', items{1}];
	[y, fs, nbits]=wavread(waveFile);
	fprintf('%d/%d: Press any key to hear %s\n', i, length(contents), contents{i}); pause
	plot((1:length(y))/fs, y);
	sound(y, fs);
	plot((1:length(y))/fs, y); grid on
end