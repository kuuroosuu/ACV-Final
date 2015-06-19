function wave2wave(wavFile1, wavFile2, fs2, nbits2, channel)
% wave2wave: audio file to audio file conversion
%
%	Usage:
%		wave2wave(wavFile1, wavFile2, fs2, nbits2, channel)
%
%	Description:
%		Conversion of wavFile1 to wavFile2 with new sample rate fs2 and new bit resolution nbits2.
%
%	Example:
%		waveFile1='主人下馬客在船.wav';
%		waveFile2='test.wav';
%		fs2=8000;
%		nbits2=8;
%		wave2wave(waveFile1, waveFile2, fs2, nbits2)
%		dos(['start ', waveFile2]);

%	Roger Jang, 2004xxxx

if nargin<1, selfdemo; return; end
if nargin<3, fs2=8000; end
if nargin<4, nbits2=8; end
if nargin<5, channel=1; end

[y, fs, nbits]=wavread(wavFile1);
channel2=size(y,2);
% 參數一樣，直接拷貝
if (fs==fs2) & (nbits==nbits2) & (channel2==channel)
	if ~strcmp(wavFile1, wavFile2)
		[parentDir, mainName, extName]=fileparts(wavFile1);
		if isempty(parentDir)			% a single file name
			absPath=which(wavFile1);
			copyfile(absPath, wavFile2, 'f');
		else					% a full file path	
			copyfile(wavFile1, wavFile2, 'f');
		end
	end
	return
end

% 參數不同，進行轉換
if channel==1		% Conversion from stereo to mono, if necessary (But not the other way around)
	y=y(:,1);
end
if fs~=fs2 | nbits~=nbits2,
	y2=resample(y, fs2, fs);		% Resample to fs2
else
	y2=y;
end

if strcmp(wavFile1, wavFile2)	% 若檔案一樣，稍等片刻，比較穩定
	pause(0.1);
end

%fprintf('Saved %s (fs=%d, bits=%d)\n', wavFile2, fs2, nbits2);
wavwrite(y2, fs2, nbits2, wavFile2);

% ====== Selfdemo
function selfdemo
wavFile1='主人下馬客在船.wav';
wavFile2='test.wav';
fs2=8000;
nbits2=8;
feval(mfilename, wavFile1, wavFile2, fs2, nbits2)