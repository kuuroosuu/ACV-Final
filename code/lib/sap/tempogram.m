function tgram=tempogram(wObj, prm)
% tempogram: Compute the tempogram from a given wave object

%	Roger Jang, 20110526

if nargin<1, selfdemo; return; end
if nargin<2, prm=[]; end
prm=prmSet(prm, wObj);

magSpecMat=vec2frameFeaMat(wObj.signal, prm);
powerSpecMat=20*log10(magSpecMat);
frameNum=size(powerSpecMat, 2);
flux=specFlux(powerSpecMat, prm);
trend=meanFilter(flux, prm.filterOrder);
novelty=flux-trend;
novelty(novelty<0)=0;
tgram=vec2frameFeaMat(novelty, prm.tgramPrm);
tgram=20*log10(tgram);

if prm.plot
	subplot(4,1,1);
	time=(1:length(wObj.signal))/wObj.fs;
	plot(time, wObj.signal);
	set(gca, 'xlim', [min(time), max(time)]);
%	xlabel('Time (sec)'); title(sprintf('Waveform of "%s"', strrep(waveFile, '_', '\_')));
	xlabel('Time (sec)'); title('Waveform');
	subplot(4,1,2);
	imagesc(powerSpecMat); axis xy
	title('Spectrogram');
	subplot(4,1,3);
	frameTime=frame2sampleIndex(1:frameNum-1, prm.frameSize, prm.overlap)/wObj.fs;
	plot(frameTime, flux, '.-');
	line(frameTime, trend, 'color', 'r');
	line(frameTime, novelty, 'color', 'm', 'linewidth', 2);
	set(gca, 'xlim', [min(frameTime), max(frameTime)]);
	xlabel('Time (sec)'); ylabel('Flux');
	subplot(4,1,4);
	imagesc(tgram); axis xy
end

% ====== Set default parameters
function prm=prmSet(prm, wObj)
% ====== Default prm
% === Prm for vec2frameFeaMat
def.frameSize=32/1000*wObj.fs;					% frameSize = 32 ms
def.overlap=def.frameSize-10/1000*wObj.fs;	% frameStep = 10 ms ===> 100 frames/sec
def.frameFcn='fftOneSide';
def.frameFcnPrm.fs=wObj.fs;
def.frameFcnPrm.nfft=4*def.frameSize;
% === Prm for spec. flux
def.useSpecNormalize=0;
def.useHalfRectifier=1;
% === Prm for obtaining novelty curve
def.filterOrder=7;
% === Prm for obtaining tempogram from novelty curve
def.tgramPrm.frameSize=2*wObj.fs/(def.frameSize-def.overlap);	% 2 sec for frame duration
def.tgramPrm.overlap=def.tgramPrm.frameSize*99/100;
def.tgramPrm.frameFcn='fftOneSide';
def.tgramPrm.frameFcnPrm.fs=wObj.fs/(def.frameSize-def.overlap);	% frame rate
def.tgramPrm.frameFcnPrm.nfft=40*def.tgramPrm.frameSize;			% zero-padding to this size.
% === For plotting
def.plot=0;
% ====== Set to the default prm if given an empty struct
if isempty(prm), prm=def; end
% ====== Set separate fields to the default values if not fully specified.
fieldName=fieldnames(def);
for i=1:length(fieldName)
	if ~isfield(prm, fieldName{i})
		prm.(fieldName{i})=def.(fieldName{i});
	end
end

% ====== Self demo
function selfdemo
waveFile='song01.wav';
wObj=waveFile2Obj(waveFile);
prm.plot=1;
flux=tempogram(wObj, prm);
