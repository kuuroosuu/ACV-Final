function specgram=wave2specgram(wObj, opt, showPlot)
% wave2specgram: Compute spectrogram from a wave object
%
%	Usage:
%		specgram=wave2specgram(wObj)
%		specgram=wave2specgram(wObj, opt)
%		specgram=wave2specgram(wObj, opt, showPlot)
%
%	Example:
%		waveFile='twinkle_twinkle_little_star.wav';
%		wObj=waveFile2obj(waveFile);
%		opt=wave2specgram('defaultOpt');
%		opt.type='log';
%		spec=wave2specgram(wObj, opt, 1);

%	Roger Jang, 20130909

if nargin<1, selfdemo; return; end
if nargin<2 && ischar(wObj) && strcmpi(wObj, 'defaultOpt')
	specgram.frameDuration=32;		% ms
	specgram.frameShift=specgram.frameDuration/2;
	specgram.type='linear';	% 'linear' (default) or 'log';
	specgram.paddedFactor=1;	% 1 for no zero padding, 2 for zero padding to have a double frame size, etc
	return
end
if nargin<3, showPlot=0; end

% === In case some fields are not defined...
if ~isfield(opt, 'type'), opt.type='linear'; end
if ~isfield(opt, 'paddedFactor'), opt.paddedFactor=1; end

wObj.signal=mean(wObj.signal, 2);	% Convert to mono
frameSize=round(wObj.fs*opt.frameDuration/1000);
overlap=round(wObj.fs.*(opt.frameDuration-opt.frameShift)/1000);

specgram.fs=wObj.fs;
specgram.frameSize=frameSize;
specgram.frameSizePadded=opt.paddedFactor*frameSize;
specgram.signal=abs(spectrogram(wObj.signal, hamming(frameSize), overlap, specgram.frameSizePadded, wObj.fs));
if strcmp(opt.type, 'log')
	specgram.signal=log(specgram.signal);
end
specgram.freqVec=linspace(0, wObj.fs/2, frameSize/2+1);
frameNum=size(specgram.signal, 2);
specgram.frameTime=frame2sampleIndex(1:frameNum, frameSize, overlap)/wObj.fs;

% === Plotting
if showPlot
	subplot(211);
	time=(0:length(wObj.signal)-1)/wObj.fs;
	plot(time, wObj.signal); set(gca, 'xlim', [min(time), max(time)]);
	xlabel('Time (sec)'); ylabel('Amplitude'); title('Waveform');
	subplot(212);
	imagesc(specgram.frameTime, specgram.freqVec, specgram.signal); axis xy
	xlabel('Time (sec)'); ylabel('Frequency (Hz)'); title('Spectrogram');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
