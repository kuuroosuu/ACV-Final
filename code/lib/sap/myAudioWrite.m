function myAudioWrite(au, audioFile)
% myAudioWrite: Save audio object to an audio file
%
%	Usage:
%		myAudioWrite(au, audioFile)
%
%	Example:
%		audioFile='what_movies_have_you_seen_recently.wav';
%		au=myAudioRead(audioFile);
%		tempFile=[tempname, '.wav'];
%		myAudioWrite(au, tempFile)
%		dos(['start ', tempFile]);
%
%	See also myAudioRead.

%	Roger Jang, 20150319

if nargin<1, selfdemo; return; end

if ~au.amplitudeNormalized
	au.signal=au.signal/(2^au.nbits/2);
end
% Old version of MATLAB which does not support audiowrite().
if verLessThan('matlab', '8.1')
	wavwrite(au.signal, au.fs, au.nbits, audioFile);
	return
end
% New version of MATLAB which supports audioread().
audiowrite(audioFile, au.signal, au.fs, 'BitsPerSample', au.nbits);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
