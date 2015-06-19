function speakerData=speakerDataRead(dirName, maxUtterNumPerSpeaker, maxSpeakerNum);
% speakerDataRead: Read wave file info for speaker identification from a given directory.
%	Usage: speakerData=speakerDataRead(dirName, maxUtterNumPerSpeaker, maxSpeakerNum);
%		dirName: root directory that contains subdirectory of each speaker, in which the wave files reside.
%		senenceNumPerSpeaker: How many sentences to read for each speaker
%		speakerData: the retrieved data structure

%	Roger Jang, 20030509, 20070517

if nargin<1, dirName='D:\dataset\tangPoem\2003-音訊處理與辨識'; end
if nargin<2, maxUtterNumPerSpeaker=inf; end
if nargin<3, maxSpeakerNum=inf; end

if (dirName(end)~='/') | (dirName(end)~='\'); dirName=[dirName, '/']; end

% ====== Collect feature from all wave files
speakerData = dir(dirName);
if length(speakerData)<2
	error(sprintf('Cannot find enough speakers/directories under "%s"!\n', dirName));
end
speakerData(1:2) = [];				% Get rid of '.' and '..'
speakerData=speakerData([speakerData.isdir]);	% Take directories only
speakerNum=min(length(speakerData), maxSpeakerNum);
speakerData=speakerData(1:speakerNum);
for i=1:speakerNum
	speakerData(i).path=[dirName, speakerData(i).name];
	waveFiles = dir([speakerData(i).path, '/*.wav']);
	speakerData(i).sentenceNum=length(waveFiles);
%	fprintf('%d/%d: Reading %d wave files recorded by %s\n', i, speakerNum, speakerData(i).sentenceNum, speakerData(i).name);
	for j=1:speakerData(i).sentenceNum
		if j>maxUtterNumPerSpeaker, break; end
		speakerData(i).sentence(j).path = [dirName, speakerData(i).name, '/', waveFiles(j).name];
%		fprintf('\t%d/%d: Processing %s...\n', j, speakerData(i).sentenceNum, speakerData(i).sentence(j).path);
		speakerData(i).sentence(j).mainName=waveFiles(j).name(1:end-4);
		items=split(speakerData(i).sentence(j).mainName, '#');
		speakerData(i).sentence(j).text=items{1};
		speakerData(i).sentence(j).speaker=speakerData(i).name;
	end
end

speakerData=rmfield(speakerData, 'isdir');	% get rid of 'isdir' field
speakerData=rmfield(speakerData, 'bytes');	% get rid of 'bytes' field