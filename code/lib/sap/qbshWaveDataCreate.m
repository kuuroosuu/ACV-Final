function waveData=qbshWaveDataCreate(waveList, qbshOpt, fileNameType)
% waveDataCreate: Create waveData from waveList for QBSH
%	Usage:
%		waveData=waveDataCreate(waveList, qbshOpt, fileNameType)
%			fileNameType=0 ===> 00028.wav
%			fileNameType=1 ===> 十個印第安人_不詳_0.wav

%	Roger Jang, 20090925

if nargin<2|isempty(qbshOpt), qbshOpt=qbshOptSet; end
if nargin<3, fileNameType=0; end

if exist(waveList)==2
	contents=textread(waveList, '%s', 'delimiter', '\n', 'whitespace', '');
	for i=1:length(contents), waveData(i).path=contents{i}; end
elseif exist(waveList)==7
	waveData=recursiveFileList(waveList, 'wav');
else
	error('Given argument is not a file or a directory!');
end
if isempty(waveData), error('Cannot read any wave files from %s!\n', waveList); end
for i=1:length(waveData)
	[parentDir, waveData(i).songName, extName]=fileparts(waveData(i).path);
	[parentDir, waveData(i).singer]=fileparts(parentDir);
end

% ====== This is for MIR-QBSH corpus only
for i=1:length(waveData)
	if fileNameType==1
		index=find(waveData(i).name=='_');	% 十個印第安人_不詳_0.wav
		waveData(i).songName=waveData(i).name(1:index(2)-1);
	end
	waveData(i).songId=waveData(i).songName;
end

% ====== Pitch tracking
%fprintf('PT over %d wav files using %s,,,\n', length(waveData), qbshOpt.ptOpt.mainFun);
for i=1:length(waveData)
	theTic=tic;
	pvFile=[waveData(i).path(1:end-3), 'pv'];
	if exist(pvFile), waveData(i).tPitch=asciiRead(pvFile); end
	% Pitch tracking
	waveData(i).cPitch=[]; waveData(i).clarity=[];
	if ~qbshOpt.usePv
		wObj=waveFile2obj(waveData(i).path);
		if isfield(qbshOpt, 'duration')
			duration=min(qbshOpt.duration, length(wObj.signal)/wObj.fs);
			wObj.signal=wObj.signal(1:floor(duration*wObj.fs));
		end
		[waveData(i).cPitch, waveData(i).clarity]=pitchTracking(waveData(i).path, qbshOpt.ptOpt);
	%	waveData(i).cPitch=pitchTrackingForcedSmooth(waveData(i).path);
	end
%	if mod(i, 100)==0
%		fprintf('%d/%d: PT of %s, time=%.4f sec\n', i, length(waveData), waveData(i).path, toc(theTic));
%	end
end

return

% ====== Check if pitch is all zero
deleteIndex=[];
for i=1:length(waveData)
	pitchSum=sum(waveData(i).pitch);
	if pitchSum==0, deleteIndex=[deleteIndex, i]; end
end
if ~isempty(deleteIndex)
	waveData(deleteIndex)=[];
	fprintf('%d wave files have all-zero pitch: %s\n', length(deleteIndex), mat2str(deleteIndex));
end