function volume=wave2volume(wObj, opt, showPlot)
% wave2volume: Wave to volume conversion
%
%	Usage:
%		volume=wave2volume(wObj, opt, showPlot)
%		wObj: wave signal
%		opt: options for volume computation
%			opt.frameSize: frame size
%			opt.overlap: overlap
%			opt.frame2volumeOpt: see the option of "frame2volume"
%
%	Example:
%		waveFile='lately2.wav';
%		wObj=waveFile2obj(waveFile);
%		opt=wave2volume('defaultOpt');
%		opt.frameSize=640; opt.overlap=480;
%		time=(1:length(wObj.signal))/wObj.fs;
%		subplot(3,1,1);
%		plot(time, wObj.signal); axis tight;
%		xlabel('Time (sec)'); ylabel('Amplitude'); title(waveFile);
%		subplot(3,1,2);
%		opt1=opt; opt1.frame2volumeOpt.method='absSum';
%		wave2volume(wObj, opt1, 1);
%		ylabel(opt1.frame2volumeOpt.method);
%		subplot(3,1,3);
%		opt2=opt; opt2.frame2volumeOpt.method='decibel';
%		wave2volume(wObj, opt2, 1);
%		ylabel(opt2.frame2volumeOpt.method);

if nargin<1, selfdemo; return; end
if ischar(wObj) && strcmp(lower(wObj), lower('defaultOpt'))
	volume.frameSize=256;
	volume.overlap=0;
	volume.frame2volumeOpt=frame2volume('defaultOpt');
	return
end
if nargin<2, opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name

frameMat=enframe(wObj.signal, opt.frameSize, opt.overlap);
volume=frame2volume(frameMat, opt.frame2volumeOpt);

if showPlot
	frameNum=size(frameMat, 2);
	frameTime=frame2sampleIndex(1:frameNum, opt.frameSize, opt.overlap)/wObj.fs;
	plot(frameTime, volume, '.-');
	ylabel(opt.frame2volumeOpt.method);
	xlabel('Time (sec)');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
