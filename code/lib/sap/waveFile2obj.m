function wObj=waveFile2obj(waveFile)
% waveFile2obj: Wave file to wave obj conversion
%
%	Usage:
%		wObj=waveFile2obj(waveFile)
%
%	Example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		wObj=waveFile2obj(waveFile);
%		disp(wObj);
%
%	See also wObj2file.

%	Roger Jang, 20100410

if nargin<1, selfdemo; return; end

if isstr(waveFile)
	[parentDir, mainName, extName]=fileparts(waveFile);
	switch lower(extName)
		case {'.wav', '.au', '.mp3'}
			if verLessThan('matlab', '8.1')
				[wObj.signal, wObj.fs, wObj.nbits]=wavread(waveFile);
			else
				info=audioinfo(waveFile);
				wObj.nbits=16;	% Default value for mp3
				if isfield(info, 'BitsPerSample'), wObj.nbits=info.BitsPerSample; end
				[wObj.signal, wObj.fs]=audioread(waveFile);
			end
	%	case '.au'
	%		[wObj.signal, wObj.fs, wObj.nbits]=auread(waveFile);
	%	case '.mp3'
	%		[wObj.signal, wObj.fs, wObj.nbits]=mp3read(waveFile);
		case '.aif'
			if verLessThan('matlab', '8.0')
				[wObj.signal, wObj.fs]=aiffread(waveFile);
			else
				[wObj.signal, wObj.fs]=audioread(waveFile);
			end
			wObj.nbits=16;
		otherwise
			error('Unknown audio file with extension=%s!', extName);
	end
end
wObj.file=waveFile;
wObj.amplitudeNormalized=1;		% The amplitude is normalized to [-1, 1] already

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
