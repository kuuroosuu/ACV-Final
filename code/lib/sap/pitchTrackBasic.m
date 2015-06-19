function [pitch, pitch0]=pitchTrackBasic(au, opt, showPlot)
% pitchTrackBasic: a simple pitch tracking method
%
%	Usage:
%		pitch=pitchTrackBasic(au)
%		pitch=pitchTrackBasic(au, opt)
%		pitch=pitchTrackBasic(au, opt, showPlot)
%
%	Description:
%		pitch=pitchTrackBaic(au, opt) return the pitch (in semitone) of an audio file by using basic frame-based PDF.
%
%	Example:
%		waveFile='soo.wav';
%		au=myAudioRead(waveFile);
%		opt=pitchTrackBasic('defaultOpt');
%		opt.frame2pitchOpt.pdf='acf';
%		opt.frame2pitchOpt.fs=au.fs;
%		showPlot=1;
%		pitch=pitchTrackBasic(au, opt, showPlot);
%
%	See also pitchTrack.

%	Category: Pitch tracking
%	Roger Jang, 20150408


if nargin<1, selfdemo; return; end
if nargin==1 && ischar(au) && strcmpi(au, 'defaultOpt')	% Set default options
	pitch.frameDuration=32;		% Duration (in ms) of a frame
	pitch.overlapDuration=0;	% Duration (in ms) of overlap
	pitch.frame2pitchOpt=frame2pitch('defaultOpt');
	pitch.useVolThreshold=1;	% 1 if the pitch is volume thresholded
	pitch.medianFilterOrder=0;	% 0 if no median filter is used to smooth pitch
	return
end
if nargin<2|isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

if isstr(au), au=myAudioRead(au); end	% If the give au is a file name
y=au.signal; fs=au.fs;
y=y-mean(y);
frameSize=round(opt.frameDuration*fs/1000);
overlap=round(opt.overlapDuration*fs/1000);
frameMat=enframe(y, frameSize, overlap);
frameNum=size(frameMat, 2);
pitch=zeros(1, frameNum);
for i=1:frameNum
	frame=frameMat(:, i);
	pitch(i)=frame2pitch(frame, opt.frame2pitchOpt);
end
volume=frame2volume(frameMat);
volumeTh=max(volume)/10;
pitch0=pitch;

if opt.medianFilterOrder>1
	pitch=medianFilter(pitch, opt.medianFilterOrder);
end
if opt.useVolThreshold
	pitch(volume<volumeTh)=nan;	% Volume-thresholded pitch
end

if showPlot
	frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/fs;
	subplot(3,1,1);
	plot((1:length(y))/fs, y); set(gca, 'xlim', [-inf inf]);
	title('Waveform');
	subplot(3,1,2);
	plot(frameTime, volume, '.-'); set(gca, 'xlim', [-inf inf]);
	line([0, length(y)/fs], volumeTh*[1, 1], 'color', 'r');
	title('Volume');
	subplot(3,1,3);
	pvFile=[au.file(1:end-3), 'pv'];
	tPitch=nan*pitch;
	if exist(pvFile, 'file')
		tPitch=asciiRead(pvFile);
		tPitch(tPitch==0)=nan;
		plot(frameTime, tPitch, 'og', frameTime, pitch0, '.-b', frameTime, pitch, '.-k');
		title('Blue: original pitch, black: volume-thresholded pitch, green: target pitch)');
	else
		plot(frameTime,	 pitch0, '.-b' , frameTime, pitch, '.-k');
		title('Blue: original pitch, black: volume-thresholded pitch)');
	end
	axis tight;
	xlabel('Time (second)'); ylabel('Semitone');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
