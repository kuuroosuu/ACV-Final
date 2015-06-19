function wav2rawFile(wavFile, rawFile)
% wav2rawFile: Convert a wav file to a raw (pcm) file

%	Roger Jang, 20081007

[y, fs, nbits]=wavReadInt(wavFile);

if size(y,2)~=1
	error('Can only handle mono wav file!');
end

if nbits==8
	rawWrite(y, rawFile, 'uint8');
elseif nbits==16
	rawWrite(y, rawFile, 'int16');
end