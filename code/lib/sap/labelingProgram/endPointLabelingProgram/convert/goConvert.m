clear all;

waveDir='../../waveFile';
waveFiles=recursiveFileList(waveDir, 'wav');
waveNum=length(waveFiles);

if waveNum>0
	fprintf('Convert %d wave files!\n', waveNum);
else
	error(sprintf('Cannot find any wave file under %s!\n', waveDir));
end

for i=1:waveNum
	source=waveFiles(i).path;
	[parentDir, mainName, extName, junk]=fileparts(waveFiles(i).path);
	epFile=[parentDir, '\endPoint.txt'];
	ep=endPointRead(epFile);
	allFileName={ep.fileName};
	index=findCellStr(upper(allFileName), upper(waveFiles(i).name));
	if length(index)~=1
		msg=sprintf('========== Something wrong with %s!', waveFiles(i).path);
		fprintf(msg);
		fid=fopen('error.txt', 'a'); fprintf(fid, '%s\n', msg); fclose(fid);
	end
	target=sprintf('%s/%s_%d_%d.wav', parentDir, mainName, ep(index).startIndex, ep(index).endIndex);
	movefile(source, target, 'f');
	fprintf('%d/%d: Move %s to %s\n', i, waveNum, source, target);
end