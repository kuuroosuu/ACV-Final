function binWave = addCueLabel2binWave(binWave, cuePosition, cueLabel)
% addCueLabel2binWave: add cue label to binary waveform
%	Usage: binWave = addCueLabel2binWave(binWave, cuePosition, cueLabel)
%       binWave:        binary waveform
%		cuePosition:    sample index of each cue point
%       cueLabel:       label of each cue point, cell array of size cue 
%           number X 1, and each element is a string; if this is empty, then
%           add cue points with default labels (e.g., Cue1, Cue2, ...) into 
%           binary waveform
% 
%   For more details about wave format, please refer to 
%   http://www.sonicspot.com/guide/wavefiles.html
%
%	Jia-Min Ren, 20100524

if nargin < 1, selfdemo, return; end
if nargin < 2, return; end
if nargin < 3, cueLabel=[]; end

cueNum = numel(cuePosition); % # of cue labels

% Cue Data Chunk Format (litter-endian format)
% Offset 	Size 	Description 	Value
% 0x00 	4 	Chunk ID 	"cue " (0x63756520), big endian
% 0x04 	4 	Chunk Data Size 	depends on Num Cue Points, little endian
% 0x08 	4 	Num Cue Points 	number of cue points in list,  little endian
% ChunkID=swapbytes(uint32(hex2dec('63756520'))); % "cue "
ChunkID = double('cue ');
% ChunkDataSize = 4 + (NumCuePoints * 24)
ChunkDataSize = dec2hex(swapbytes((uint32(4+cueNum*24))));
% Note that in the case of the bits of variable is less than 8, it must add
% zero bit before the first bit
zeroPadding = repmat('0', 1, 8-numel(ChunkDataSize));
ChunkDataSize = [zeroPadding, ChunkDataSize];
ChunkDataSize = [hex2dec(ChunkDataSize(1:2)), hex2dec(ChunkDataSize(3:4)),...
    hex2dec(ChunkDataSize(5:6)), hex2dec(ChunkDataSize(7:8))];
NumCuePoints=dec2hex(swapbytes((uint32(cueNum))));
% Adding zero bits if necessary
zeroPadding = repmat('0', 1, 8-numel(NumCuePoints));
NumCuePoints = [zeroPadding, NumCuePoints];
NumCuePoints = [hex2dec(NumCuePoints(1:2)), hex2dec(NumCuePoints(3:4)),...
    hex2dec(NumCuePoints(5:6)), hex2dec(NumCuePoints(7:8))];
% write to binary wave content
binWave = [binWave; ChunkID'; ChunkDataSize'; NumCuePoints'];

fileSize = 12;

% List of Cue Points
% A list of cue points is simply a set of consecutive cue point descriptions 
% that follow the format described below.
% Offset 	Size 	Description 	Value
% 0x00 	4 	ID 	unique identification value, little endian
% 0x04 	4 	Position 	play order position, little endian
% 0x08 	4 	Data Chunk ID 	RIFF ID of corresponding data chunk, big endian
% 0x0c 	4 	Chunk Start 	Byte Offset of Data Chunk *
% 0x10 	4 	Block Start 	Byte Offset to sample of First Channel
% 0x14 	4 	Sample Offset 	Byte Offset to sample byte of First Channel, little endian
CuePoints=[];
for cueIdx=1:cueNum,
    ID=dec2hex(swapbytes((uint32(cueIdx))));
    % Adding zero bits if necessary
    zeroPadding = repmat('0', 1, 8-numel(ID));
    ID = [zeroPadding, ID];
    ID=[hex2dec(ID(1:2)), hex2dec(ID(3:4)), hex2dec(ID(5:6)), hex2dec(ID(7:8))];    
    Position=dec2hex(swapbytes((uint32(cuePosition(cueIdx)))));
    % Adding zero bits if necessary
    zeroPadding = repmat('0', 1, 8-numel(Position));
    Position = [zeroPadding, Position];
    Position=[hex2dec(Position(1:2)), hex2dec(Position(3:4)), hex2dec(Position(5:6)), hex2dec(Position(7:8))];        
    DataChunkID = double('data');
    ChunkStart = zeros(1,4);
    BlockStart = zeros(1,4);
    SampleOffset = Position;
    CuePoints=[CuePoints, ID, Position, DataChunkID, ChunkStart, BlockStart, SampleOffset];
end
binWave = [binWave; CuePoints'];

fileSize = fileSize + 24*cueNum;

% Associated Data List Chunk - "LIST"
% Offset 	Size 	Description 	Value
% 0x00 	4 	Chunk ID 	"LIST" (0x4C495354),  big endian
% 0x04 	4 	Chunk Data Size 	depends on contained text, little endian
% 0x08 	4 	Type ID 	"adtl" (0x6164746C), big endian
% 0x0c 	List of Text Labels and Names
% ChunkID=swapbytes(uint32(hex2dec('4c495354'))); % list
ChunkID = double('LIST');
ChunkDataSize = 4 + (28*cueNum);
% if # of chars in a cuelabel is odd, then add data size to # of chars + 1;
% else (# of chars is even) then add data size to # of chars + 2
sum=0;
for idx = 1 : numel(cueLabel),
    label=cueLabel{idx};
    if ~isempty(label),
        sum=sum+12;
        if mod(numel(label),2)==0,
            sum=sum+numel(label)+2;
        else
            sum=sum+numel(label)+1;
        end
    end
end

ChunkDataSize = dec2hex(swapbytes((uint32(ChunkDataSize + sum))));
% Adding zero bits if necessary
zeroPadding = repmat('0', 1, 8-numel(ChunkDataSize));
ChunkDataSize = [zeroPadding, ChunkDataSize];
ChunkDataSize = [hex2dec(ChunkDataSize(1:2)), hex2dec(ChunkDataSize(3:4)),...
    hex2dec(ChunkDataSize(5:6)), hex2dec(ChunkDataSize(7:8))];
% TypeID= swapbytes(uint32(hex2dec('6164746C'))); %double('adtl');
TypeID = double('adtl');
binWave = [binWave; ChunkID'; ChunkDataSize'; TypeID'];

fileSize = fileSize+12;

% Labeled Text Chunk - "ltxt"
% Offset 	Size 	Description 	Value
% 0x00 	4 	Chunk ID 	"ltxt" (0x6C747874),                big endian
% 0x04 	4 	Chunk Data Size 	depends on contained text,  little endian
% Note that Chunk Data Size is always (0x14000000) in the sample
% 0x08 	4 	Cue Point ID 	0 - 0xFFFFFFFF,                 little endian
% 0x0c 	4 	Sample Length 	0 - 0xFFFFFFFF,                 little endian
% 0x10 	4 	Purpose ID 	0 - 0xFFFFFFFF,(0x72676E20),"rgn ", little endian
% 0x12 	2 	Country 	0 - 0xFFFF,                         little endian
% 0x14 	2 	Language 	0 - 0xFFFF,                         little endian
% 0x16 	2 	Dialect 	0 - 0xFFFF,                         little endian
% 0x18 	2 	Code Page 	0 - 0xFFFF,                         little endian
ltxt=[];
for cueIdx = 1 : cueNum,
    ChunkID = double('ltxt');
    ChunkDataSize = [hex2dec('14'), hex2dec('00'), hex2dec('00'), hex2dec('00')];
    ID=dec2hex(swapbytes((uint32(cueIdx))));
    % Adding zero bits if necessary
    zeroPadding = repmat('0', 1, 8-numel(ID));
    ID = [zeroPadding, ID];
    ID=[hex2dec(ID(1:2)), hex2dec(ID(3:4)), hex2dec(ID(5:6)), hex2dec(ID(7:8))];
    SampleLength = zeros(1,4);
    PurposeID = double('rgn ');
    CountryLanguage = zeros(1,4);
    DialectCodePage = zeros(1,4);
    ltxt = [ltxt, ChunkID, ChunkDataSize, ID, SampleLength, PurposeID, CountryLanguage, DialectCodePage];
end

binWave = [binWave; ltxt'];

fileSize=fileSize+28*cueNum;

% Label Chunk - "labl"
% Offset 	Size 	Description 	Value
% 0x00 	4 	Chunk ID 	"labl" (0x6C61626C),                big endian
% 0x04 	4 	Chunk Data Size 	depends on contained text,  little endian
% 0x08 	4 	Cue Point ID 	0 - 0xFFFFFFFF,                 little endian
% 0x0c 	Text,                                               big endian
% Note that 
% The text is a null terminated string of characters. If the number of 
% characters in the string is not even, padding must be appended to the 
% string. The appended padding is not considered in the label chunk's chunk
% size field. 

% ChunkID = swapbytes(uint32(hex2dec('6c61626c')));
ChunkID = double('labl');

labl=[];
% % label
% for cueIdx = 1 : numel(cueLabel),
%     label = cueLabel{cueIdx};
%     % ignore empty cue label
%     if isempty(label), continue; end,
%     ChunkDataSize = uint32(5 + numel(label));
%     CuePointID = uint32(cueIdx);
% 
%     text = [];
%     charNum = numel(label);
%     itr = floor(charNum / 4);
%     for idx=1:itr,
%         tmp = dec2hex(double(label((idx-1)*4+1:idx*4)))';
%         tmp = tmp(:)';
%         text(idx) = swapbytes(uint32(hex2dec(tmp)));
%     end
%     
%     remainingCharNum = mod(charNum,4);
%     switch remainingCharNum,
%         case 1,
%             tmp = dec2hex(double(label(end)))';
%             tmp = [tmp(:)' '000000'];
%             text(end+1) = swapbytes(uint32(hex2dec(tmp)));
%         case 2,
%             tmp = dec2hex(double(label(end-1:end)))';
%             tmp = [tmp(:)' '0000'];
%             text(end+1) = swapbytes(uint32(hex2dec(tmp)));
%         case 3,
%             tmp = dec2hex(double(label(end-2:end)))';
%             tmp = [tmp(:)' '00'];
%             text(end+1) = swapbytes(uint32(hex2dec(tmp)));            
%     end
% 
%     labl = [labl; ChunkID; ChunkDataSize; CuePointID; text'];
%     
%     fileSize = fileSize + 12; % ChunkID, ChunkDataSize, CuePointID
%     fileSize = fileSize + 4*numel(text);    % label data
%     
% end

% label
for cueIdx = 1 : numel(cuePosition),
    if isempty(cueLabel),
        label = ['Cue' num2str(cueIdx)];
    else
        label = cueLabel{cueIdx};
    end
    % ignore empty cue label
    if isempty(label), continue; end,
    ChunkDataSize = dec2hex(swapbytes((uint32(5 + numel(label)))));
    % Adding zero bits if necessary
    zeroPadding = repmat('0', 1, 8-numel(ChunkDataSize));
    ChunkDataSize = [zeroPadding, ChunkDataSize];
    ChunkDataSize = [hex2dec(ChunkDataSize(1:2)), hex2dec(ChunkDataSize(3:4)),...
        hex2dec(ChunkDataSize(5:6)), hex2dec(ChunkDataSize(7:8))];    
    ID=dec2hex(swapbytes((uint32(cueIdx))));
    % Adding zero bits if necessary
    zeroPadding = repmat('0', 1, 8-numel(ID));
    ID = [zeroPadding, ID];
    ID=[hex2dec(ID(1:2)), hex2dec(ID(3:4)), hex2dec(ID(5:6)), hex2dec(ID(7:8))];

    if mod(numel(label),2)==0, % even characters, then add 00
        text = [double(label),0,0];
    else % odd characters, add 0
        text = [double(label),0];
    end

    labl = [labl, ChunkID, ChunkDataSize, ID, text];
    
    fileSize = fileSize + 12; % ChunkID, ChunkDataSize, CuePointID
    fileSize = fileSize + numel(text);    % label data
    
end
binWave = [binWave; labl'];

% process file size
fileSize = fileSize + binWave(5)+2^8*binWave(6)+2^16*binWave(7)+2^24*binWave(8);
fileSize=dec2hex(swapbytes((uint32(fileSize))));
% Adding zero bits if necessary
zeroPadding = repmat('0', 1, 8-numel(fileSize));
fileSize = [zeroPadding, fileSize];
fileSize=[hex2dec(fileSize(1:2)), hex2dec(fileSize(3:4)), hex2dec(fileSize(5:6)), hex2dec(fileSize(7:8))];
binWave(5:8) = fileSize;


function selfdemo
waveFile = 'test.wav';
fid = fopen(waveFile, 'rb');
binWave = fread(fid); fclose(fid);
cuePosition=[5000,10000,12000];
cueLabel={'U12','U23','U12'};
binWave = feval(mfilename, binWave, cuePosition, cueLabel);
waveFile = 'testCue.wav';
fid = fopen(waveFile,'wb');
fwrite(fid,binWave); fclose(fid);
dos('start testCue.wav');