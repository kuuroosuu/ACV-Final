clear all;
pdOpt=pdOptSet;

fprintf('Collecting wave data and features from "%s"...\n', pdOpt.waveDir);
waveData=recursiveFileList(pdOpt.waveDir, 'wav');	% Collect all wave files

totalCount=[0 0 0 0]';
for i=1:length(waveData)
	fprintf('%d/%d ==> %s\n', i, length(waveData), waveData(i).name);
	waveFile=waveData(i).path;
	pvFile=strrep(waveFile, '.wav', '.pv');
	if ~exist(pvFile)
		fprintf('\tNo .pv file found for %s\n', waveFile);
		continue;
	end
	pv=asciiRead(pvFile);
	classVec=pv>0;
	mat=buffer2(classVec, 2, 1)';
	[sortedRow, theCount]=rowCount(mat);
	totalCount=totalCount+theCount;
end

transCount=reshape(totalCount, 2, 2)';
transProb=transCount./(sum(transCount, 2)*ones(1,2));
save transProb transProb