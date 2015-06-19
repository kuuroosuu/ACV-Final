function au=myAudioRead(audioFile)
% myAudioRead: Read an audio file to create an audio object
%
%	Usage:
%		au=myAudioRead(audioFile)
%
%	Example:
%		audioFile='youAtLeast.mp3';
%		au=myAudioRead(audioFile);
%		disp(au);
%
%	See also audio2file.

%	Roger Jang, 20141214

if nargin<1, selfdemo; return; end

[parentDir, mainName, extName]=fileparts(audioFile);
% Old version of MATLAB which does not support audioread().
if verLessThan('matlab', '8.1')
	switch lower(extName)
		case '.wav'
			[au.signal, au.fs, au.nbits]=wavread(audioFile);
		case '.au'
			[au.signal, au.fs, au.nbits]=auread(audioFile);
		case '.mp3'
			[au.signal, au.fs, au.nbits]=mp3read(audioFile);
		case '.aif'
			[au.signal, au.fs]=aiffread(audioFile);
			au.fs=16;		% Sure?
	end
	au.file=audioFile;
	au.amplitudeNormalized=1;		% The amplitude is normalized to [-1, 1] already
	return
end

% New version of MATLAB which supports audioread().
[au.signal, au.fs]=audioread(audioFile);
info=audioinfo(audioFile);
if isfield(info, 'BitsPerSample')
	au.nbits=info.BitsPerSample;
else
	fprintf('Cannot identify BitsPerSample of %s, assuming it is 16 bits\n', audioFile)
	au.nbits=16;
end

au.file=audioFile;
au.amplitudeNormalized=1;		% The amplitude is normalized to [-1, 1] already

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
