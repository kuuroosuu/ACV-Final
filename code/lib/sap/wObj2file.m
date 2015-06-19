function wObj2file(wObj, waveFile)
% wObj2file: Save wave object to a wave file
%
%	Usage:
%		wObj2file(wObj, waveFile)
%
%	Example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		wObj=waveFile2obj(waveFile);
%		tempWaveFile=[tempname, '.wav'];
%		wObj2file(wObj, tempWaveFile)
%		dos(['start ', tempWaveFile]);
%
%	See also waveFile2obj.

%	Roger Jang, 20120607

if nargin<1, selfdemo; return; end

if ~wObj.amplitudeNormalized
	wObj.signal=wObj.signal/(2^wObj.nbits/2);
end
%wavwrite(wObj.signal, wObj.fs, wObj.nbits, waveFile);
audiowrite(waveFile, wObj.signal, wObj.fs);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
