function pitch=pitchTrackingSimple(waveFile, opt, showPlot)
% pitchTrackingSimple: a simple pitch tracking method
%
%	Usage:
%		pitch=pitchTrackingSimple(waveFile, opt, plotOpt)
%
%	Example:
%		waveFile='soo.wav';
%		opt.frameDuration=32;	% Duration (in ms) of a frame
%		opt.overlapDuration=0;	% Duration (in ms) of overlap
%		opt.maxShiftDuration=32;
%		opt.pdf='frame2acf';
%		opt.pdfMethod=1;
%		showPlot=1;
%		pitch=pitchTrackingSimple(waveFile, opt, showPlot);

if nargin<1, selfdemo; return; end
if nargin==1 && ischar(waveFile) && strcmpi(waveFile, 'defaultOpt')	% Set default options
	pitch.frame2pdfOpt=frame2pdf('defaultOpt');
	pitch.frameDuration=32;		% Duration (in ms) of a frame
	pitch.overlapDuration=0;	% Duration (in ms) of overlap
	pitch.maxShiftDuration=32;
	return
end
if nargin<2|isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

% ====== Sync options

[y, fs]=wavread(waveFile);
y=y-mean(y);
frameSize=round(opt.frameDuration*fs/1000);
overlap=round(opt.overlapDuration*fs/1000);
opt.frame2pdfOpt.maxShift=round(opt.maxShiftDuration*fs/1000);
maxFreq=1000;
minFreq=40;
n1=round(fs/maxFreq);	% pdf(1:n1) will not be used
n2=round(fs/minFreq);	% pdf(n2:end) will not be used
frameMat=enframe(y, frameSize, overlap);
frameNum=size(frameMat, 2);
volume=frame2volume(frameMat);
volumeTh=max(volume)/10;
pitch=0*volume;
for i=1:frameNum
%	fprintf('%d/%d\n', i, frameNum);
	frame=frameMat(:, i);
	pdf=frame2pdf(frame, opt.frame2pdfOpt);
	pdf(1:n1)=-inf;
	pdf(n2:end)=-inf;
	[maxValue, maxIndex]=max(pdf);
	freq=fs/(maxIndex-1);
	pitch(i)=freq2pitch(freq);
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
	pitch2=pitch; pitch2(volume<volumeTh)=nan;	% Volume-thresholded pitch
	plot(frameTime, pitch, frameTime, pitch2, '.-r'); axis tight;
	xlabel('Time (second)'); ylabel('Semitone');
	title('Original pitch (blue) and volume-thresholded pitch (red)');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
