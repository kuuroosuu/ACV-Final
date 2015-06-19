function songDb=songDbCreate(midiOrWaveList, ptOpt)
%songDbCreate: Create song database from a give midi/wave list, or from a midi directory
%	Usage: songDb=songDbCreate(midiOrWaveList)
%		songDb=songDbCreate(midiDir)

if nargin<1, midiOrWaveList='midi0048.list'; end

if exist(midiOrWaveList)==7	% This is in fact a directory
	midiData=recursiveFileList(midiOrWaveList, 'mid');
	waveData=recursiveFileList(midiOrWaveList, 'wav');
	if isempty(waveData)
		allData=midiData;
	elseif isempty(midiData)
		allData=waveData;
	else
		allData=[midiData, waveData];
	end
	contents={allData.path};
elseif exist(midiOrWaveList)==2		% A list of midi/wave files
	contents=textread(midiOrWaveList, '%s', 'delimiter', '\n', 'whitespace', '');
else
	error('The given file/folder does not exist!');
end	

fileNum=length(contents);

for i=1:fileNum
	fileName=contents{i};
	fprintf('%d/%d: Adding "%s" to song database...\n', i, fileNum, fileName);
	songDb(i).path=fileName;
	[parentDir, songDb(i).songName, songDb(i).fileType]=fileparts(songDb(i).path);
	switch lower(songDb(i).fileType)
		case '.mid'
			songDb(i).note=midiFile2note(songDb(i).path);
		case '.wav';
			songDb(i).pv=pitchTracking(songDb(i).path, ptOpt);
		otherwise
			fprintf(sprintf('Unknown file type: %s\n', songDb(i).fileType));
	end
end

% ====== This is for MIR-QBSH corpus only
for i=1:length(songDb)
	songDb(i).id=songDb(i).songName;
end

%for i=1:length(songDb)
%	if strcmp(songDb(i).fileType, '.mid')
%		[songDb(i).pv, songDb(i).noteStartIndex]=note2pv(songDb(i).note, (ptOpt.frameSize-ptOpt.overlap)/ptOpt.fs, inf, 0);
%	end
%end
