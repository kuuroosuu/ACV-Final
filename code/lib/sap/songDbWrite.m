function songDbWrite(songData, txtFile, dbFile)

if nargin<2, txtFile='test.txt'; end
if nargin<3, dbFile='test.db'; end

% ====== 寫入 txt 檔案
songNum=length(songData);
fid=fopen(txtFile, 'w');
for i=1:songNum
	fprintf(fid, 'SSN=%d\t', songData(i).SSN);
	fprintf(fid, 'songName=%s\t', songData(i).songName);
	fprintf(fid, 'language=%s\t', songData(i).language);
	fprintf(fid, 'singer=%s\t', songData(i).singer);
	fprintf(fid, 'price=%d\t', songData(i).price);
	fprintf(fid, 'calorie=%g\t', songData(i).calorie);
	fprintf(fid, 'start=%d\t', songData(i).start);
	fprintf(fid, 'size=%d\t', songData(i).size);
	fprintf(fid, 'refrain=');
	for j=1:length(songData(i).refrain)
		if j==1
			fprintf(fid, '%d', songData(i).refrain(j));
		else
			fprintf(fid, '+%d', songData(i).refrain(j));
		end
	end
	fprintf(fid, '\n');
end
fclose(fid);

% ====== 寫入 db 檔案
fid=fopen(dbFile, 'wb');
for i=1:songNum
	fwrite(fid, songData(i).track, 'integer*2');
	fwrite(fid, 255, 'integer*2');
end
fclose(fid);