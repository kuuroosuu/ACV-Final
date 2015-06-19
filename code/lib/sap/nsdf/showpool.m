function showpool
	load allscore.mat allscore
	waveDir='D:\LAB\wave\';
	waveDir=strrep(waveDir, '\', '/');
	allWaves=recursiveFileList(waveDir, 'wav');
	waveNum=length(allWaves);
	[sortscore, index] = sort(allscore);
	len = length(allscore);
	
	i=1;
	while sortscore(i)<85
		fprintf('%s...score=%g\n',allWaves(i).path,sortscore(i));
		show2pitch(allWaves(i).path, 1);
		fprintf('hit any key to next..');
		pause
		i=i+1;
	end
		
	