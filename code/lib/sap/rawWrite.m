function rawWrite(matrix, rawFile, format)
% rawWrite: Write data to a binary .raw file
%	Usage: rawWrite(matrix, rawFile)

if nargin<3, format='uint8'; end

fid=fopen(rawFile, 'w');
fwrite(fid, matrix, format);
fclose(fid);
%fprintf('Written raw data to "%s".\n', rawFile);