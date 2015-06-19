waveDir='D:\LAB\wave\';
waveDir=strrep(waveDir, '\', '/');
allWaves=recursiveFileList(waveDir, 'wav');
waveNum=length(allWaves);
fprintf('從 "%s" 讀取了 %d 個 wav 檔案\n', waveDir, waveNum);

totalTime=0;
for i=1:2
	tic;
	waveFile=allWaves(i).path;
	fprintf('%d/%d: Check the pitch of %s...\n', i, waveNum, waveFile);
	pitchLabel(waveFile);
	t=toc;
	fprintf('used time=%g\n',t);
	totalTime=totalTime+toc;
	%fprintf('Hit any key to check next wav file...\n'); pause
	%savePitch;
	%close all
end
fprintf('totalTime=%g   averageTime=%g\n',totalTime,totalTime/waveNum);