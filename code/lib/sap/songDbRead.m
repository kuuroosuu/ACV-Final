function songDb=songDbRead(dbName, encoded)
%songDbRead: Read song database files into a songDb structure
%	Usage:
%		songDb=songDbRead(dbName, encoded)
%		dbName: name of the database, which are spread in two files, dbName.txt and dbName.db
%		encoded: 0 if the database if not encoded; 1 if database is encoded
%		songDb: a structure contains fields specified in dbName.txt file
%			In particular, songDb(i).track is a vector containing (pitch, duration) pairs read from dbFile.
%			pitch is in semitone, duration is in 1/64 seconds
%
%	Example:
%		songDb=songDbRead('childSong');
%		structDispInHtml(songDb);

%	Roger Jang, 20051202, 20080526

if nargin<0, selfdemo; return; end
if nargin<1, dbName='childSong'; end
if nargin<2, encoded=0; end

% ====== 從 txtFile 讀入基本的 songDb
txtFile=[dbName, '.txt'];
dbFile=[dbName, '.db'];
songDb=tableRead(txtFile, 3);
songNum=length(songDb);
% ====== 將某幾個欄位變成數字，方便處理
for i=1:songNum
%	songDb(i).SSN=eval(songDb(i).SSN);
%	songDb(i).price=eval(songDb(i).price);
%	songDb(i).calorie=eval(songDb(i).calorie);
	songDb(i).start=eval(songDb(i).start);
	songDb(i).size=eval(songDb(i).size);
%	songDb(i).refrain=eval(['[', strrep(songDb(i).refrain, '+', ','), ']']);
end

% ====== 從 dbFile 讀入音高及因長
dbFid=fopen(dbFile,'r'); dbContents=fread(dbFid, inf, 'uint8'); fclose(dbFid);
if encoded
	dbContents=encodeXor(dbContents, 'Karami3');
end
for i=1:songNum
	startPos=songDb(i).start;
	len=songDb(i).size;
	% ====== 以下是使用 file stream 的方式來讀取
%	fseek(dbFid, startPos, 'bof');				%start position
%	songDb(i).track = fread(dbFid, len, 'uint8');	% data length
%	in = fread(dbFid, 1, 'uint8');	% check 0xFF
%	if in~=255
%		error(sprintf('[%s] data error!\n', dbFile));
%	end
	% ====== 以下是一次讀入檔案
	songDb(i).track=dbContents(startPos+1:startPos+len);
	if dbContents(startPos+len+1)~=255		% check 0xFF
		error(sprintf('"%s" data error!\n', dbFile));
	end
end

% ====== Cast to uint8 for saving memory
for i = 1:length(songDb)
	songDb(i).track = uint8(songDb(i).track);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
