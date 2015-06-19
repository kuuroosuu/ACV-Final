function wavePitchPlay(wObj, pObj, outputWaveFile)
% wavePitchPlay: Playback of wave and pitch in different channels simutaneously
%
%	Usage:
%		wavePitchPlay(wObj, pObj)
%		wavePitchPlay(wObj, pObj, outputWaveFile)
%
%	Example:
%		wObj=waveFile2obj('10LittleIndians.wav');
%		pObj.signal=asciiRead('10LittleIndians.pv');
%		pObj.frameRate=8000/256;
%		plot((1:length(wObj.signal))/wObj.fs, wObj.signal);
%		wavePitchPlay(wObj, pObj);
%
%	See also wavePitchPlayButton.

%	Category: Playback
%	Roger Jang, 20121213

if nargin<1, selfdemo; return; end

pitchLen=length(pObj.signal);
duration=ones(1, pitchLen)/pObj.frameRate;
note=[pObj.signal(:)'; duration]; note=note(:);
pitchWave=note2wave(note, 1, wObj.fs, 2);
origWave=mean(wObj.signal,2);
len=min(length(pitchWave), length(origWave)); pitchWave=pitchWave(1:len); origWave=origWave(1:len);
y=[origWave, pitchWave];
if nargin<3
	sound(y, wObj.fs);
else
%	fprintf('Saving output wave file to <a href="matlab: dos([''start '', ''%s'']);">%s</a>...\n', outputWaveFile, outputWaveFile);
	wavwrite(y, wObj.fs, 16, outputWaveFile);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
