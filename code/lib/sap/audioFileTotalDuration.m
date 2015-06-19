function [totalDuration, eachDuration]=audioFileTotalDuration(audioDir, extName)
%audioFileTotalDuration: Return the total duration of audio files in a directory
%
%	Usage:
%		[totalDuration, eachDuration]=audioFileTotalDuration(audioDir, extName)
%
%	Example:
%		audioDir='d:/users/jang/matlab/toolbox/sap';
%		extName='wav';
%		[totalDuration, eachDuration]=audioFileTotalDuration(audioDir, extName);
%		fprintf('Total duration of "%s" file in "%s" = %g sec\n', extName, audioDir, totalDuration);
%
%	See also waveFile2obj.

%	Category: String processing
%	Roger Jang, 20140801

if nargin<1, selfdemo; return; end
if nargin<2, extName='mp3'; end

audioData=recursiveFileList(audioDir, extName);
for i=1:length(audioData)
	[y, fs]=audioread(audioData(i).path);
	audioData(i).duration=size(y,1)/fs;
end
eachDuration=[audioData.duration];
totalDuration=sum(eachDuration);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);