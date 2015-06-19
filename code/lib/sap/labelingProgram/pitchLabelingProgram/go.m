% Label pitch of wave files in a given directory

addpath d:/users/jang/matlab/toolbox/utility -end
addpath d:/users/jang/matlab/toolbox/sap -end
close all; clear all;

% Directory of the wave files
auDir='D:\users\jang\books\audioSignalProcessing\programmingContest\pitchTracking\exampleProgram\waveFile\rogerJang';

auDir=strrep(auDir, '\', '/');
auSet=recursiveFileList(auDir, 'wav');
auNum=length(auSet);
fprintf('Read %d wave files from "%S"\n', auNum, auDir);

for i=1:auNum
	fprintf('%d/%d: Check the pitch of %s...\n', i, auNum, auSet(i).path);
	pitchLabel(auSet(i).path);
	fprintf('\tHit any key to check next wav file...\n'); pause
	pitchSave;
	close all
end