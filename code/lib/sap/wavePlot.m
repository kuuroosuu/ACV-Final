function wavePlot(wObj)
% wavePlot: Plot the waveform in a wav file
%
%	Usage:
%		wavePlot(waveFile)
%
%	Example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		wavePlot(waveFile);

%	Roger Jang, 20100918

if nargin<1, selfdemo; return; end

if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name

plot((1:length(wObj.signal))/wObj.fs, wObj.signal, '.-');
grid on; axis([-inf inf -1 1]);
title(strPurify4label(sprintf('waveFile=%s, fs=%d, nbits=%d', wObj.file, wObj.fs, wObj.nbits)));
xlabel('Time (sec)');
ylabel('Amplitude');
wavePlayButton(wObj);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
