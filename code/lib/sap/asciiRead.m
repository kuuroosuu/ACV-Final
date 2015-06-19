function matrix=asciiRead(fileName)
% asciiRead: Read an ascii file into a matrix
%	Usage: matrix=asciiRead(fileName)

%	Roger Jang, 20041021

if nargin<1, selfdemo; return; end

fid=fopen(fileName, 'rb');
content=fread(fid, inf, 'char');
fclose(fid);

file=tempname;
%copyfile(fileName, file);	% problematic in version 6.5!
fid=fopen(file, 'wb');
fwrite(fid, content, 'char');
fclose(fid);
load(file);
[junk, b, junk]=fileparts(file);
eval(['matrix=', b, ';']);
delete(file);

% ====== Self demo
function selfdemo
x=magic(5);
fileName='我的檔案.txt';
asciiWrite(x, fileName);
matrix=asciiRead(fileName);
fprintf('從 "%s" 讀出的內容：\n', fileName);
disp(matrix);