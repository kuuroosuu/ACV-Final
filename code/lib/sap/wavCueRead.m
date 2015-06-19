function cueTime=wavCueRead(waveFile, showPlot)
% wavCueRead: get sample indices of cue points from a given binary waveform 
%
%	Usage:
%		cueTime=wavCueRead(binWave, showPlot)
%
%	Description:
%		For more details about wave format, please refer to "http://www.sonicspot.com/guide/wavefiles.htm".
%		How to insert cue list into a wav file:
%			1. 開啟coolEdit
%			2. 利用Mouse點選欲標記的wave波型位置
%			3. 按快速鍵 F8 標記 Cue
%			4. 重複  2~3步驟   直到標記完畢
%			5. 存檔
%
%		注意：標記Cue的順序  可以不用依wave的時序
%			例如  cue1    cue3   cue2
%			wavCueRead 程式最後會做sort
%
%	Example:
%		waveFile='tapping.wav';
%		cueLabel=wavCueRead(waveFile, 1);
%		system(['start ', waveFile]);

%	Jia-Min Ren, 20090527, Roger Jang, 20150127

if nargin<1, selfdemo, return; end
if nargin<2, showPlot=0; end

% Read the wav file
fid = fopen(waveFile, 'rb');
binWave = fread(fid);
fclose(fid);

% 關鍵字'cue '在binary wave中的index...
cuePos = findstr(char(binWave'), 'cue ');
% binary格式的wave檔，一定要有關鍵字'cue '
% 若沒有，則代表waveform沒有cue points
if isempty(cuePos),
    cueLabel=[]; return;
end

%% 找出waveform裡有多少個cue points...
% Cue Data Chunk Format (litter-endian format)
% Offset 	Size 	Description 	Value
% 0x00 	4 	Chunk ID 	"cue " (0x63756520)
% 0x04 	4 	Chunk Data Size 	depends on Num Cue Points
% 0x08 	4 	Num Cue Points 	number of cue points in list
cuePos = cuePos(1); % 第一個才是真的index...(若有多個，可能是cue label名稱)
binCueNum = binWave(cuePos+8:cuePos+11);
% store in litter-endian format, so it needs read from left to right...
cueNum = binCueNum(1)+2^8*binCueNum(2)+2^16*binCueNum(3)+2^24*binCueNum(4);

%% 找出每個cue point所對應的sample index
% List of Cue Points
% A list of cue points is simply a set of consecutive cue point descriptions 
% that follow the format described below.
% Offset 	Size 	Description 	Value
% 0x00 	4 	ID 	unique identification value
% 0x04 	4 	Position 	play order position
% 0x08 	4 	Data Chunk ID 	RIFF ID of corresponding data chunk
% 0x0c 	4 	Chunk Start 	Byte Offset of Data Chunk *
% 0x10 	4 	Block Start 	Byte Offset to sample of First Channel
% 0x14 	4 	Sample Offset 	Byte Offset to sample byte of First Channel
currCuePos = cuePos + 12;  % the position of an element in the list of cue Points
cueLabel=zeros(1,cueNum);  % output data...
for cueIdx = 1 : cueNum,
    % 4 bytes, Position 	play order position
    binCuePos = binWave(currCuePos+4:currCuePos+7);
    % transform to sample index (little-endian format)
    cueLabel(cueIdx) = binCuePos(1)+2^8*binCuePos(2)+2^16*binCuePos(3)+2^24*binCuePos(4);
    % go to the position of the next cue point
    % each list of cue point occpy 24 bytes
    currCuePos = currCuePos + 4*6;
end
% embedded labels may be in any orders, so need to sort here
cueLabel=sort(cueLabel);
aObj=myAudioRead(waveFile);
y=aObj.signal; fs=aObj.fs;
cueTime=cueLabel/fs;

if showPlot
	channelNum=size(y,2);
	time=(1:size(y,1))/fs;
	if channelNum==1
		plot(time, y);
		for i=1:length(cueTime), line(cueTime(i)*[1 1], [-1 1], 'color', 'r'); end
	else
		subplot(211); plot(time, y(:,1));
		for i=1:length(cueTime), line(cueTime(i)*[1 1], [-1 1], 'color', 'r'); end
		subplot(212); plot(time, y(:,2));
		for i=1:length(cueTime), line(cueTime(i)*[1 1], [-1 1], 'color', 'r'); end
	end
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
