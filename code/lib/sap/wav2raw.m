function wav2raw(wavFile, rawFile)
% WAV2RAW Convert a wav file to a raw file

if nargin~=2,
	error('Need to input arguments');
end

[y, fs, nbits]=wavread(wavFile);
y=y*(2^nbits/2)+2^nbits/2;	% Convert to unsigned integer
if nbits==16,
	y=y/256;	% Convert to 8-bit resolution
end

fid=fopen(rawFile, 'wb');
fwrite(fid, y);		% Save to a binary file with a byte per sample
fclose(fid);
