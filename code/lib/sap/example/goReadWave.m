waveFile='sample.wav';
fid=fopen(waveFile, 'r');

n=4; x=fread(fid, 4);
fprintf('offset = %d, length = %d, content = %s\n', ftell(fid)-n, n, char(x'));

n=4; x=fread(fid, 4);
fprintf('offset = %d, length = %d, content = %d\n', ftell(fid)-n, n, x);

fclose(fid);