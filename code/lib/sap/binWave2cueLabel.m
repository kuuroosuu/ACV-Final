function cueLabel = binWave2CueLabel(binWave)
% binWave2CueLabel: get sample indices of cue points from a given binary
% waveform 
%	Usage: cueLabel = binWave2CueLabel(binWave)
%       binWave:    binary waveform (it must contain some cue points)
%		cueLabel:   sample index of each cue point
% 
%   For more details about wave format, please refer to 
%   http://www.sonicspot.com/guide/wavefiles.htm
%
%	Jia-Min Ren, 20090527

if nargin < 1, selfdemo, return; end,

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

function selfdemo
waveFile = 'testSong.wav';
fid = fopen(waveFile, 'rb');
binWave = fread(fid);
fclose(fid);
cueLabel = feval(mfilename, binWave);