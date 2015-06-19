function flux=wave2specFlux(wObj, prm)
% specFlux: Compute the flux of a given wave object

%	Roger Jang, 20110526

if nargin<1, selfdemo; return; end
if nargin<2, prm=[]; end
prm=prmSet(prm, wObj);

magSpecMat=vec2frameFeaMat(wObj.signal, prm);
powerSpecMat=20*log10(magSpecMat);
frameNum=size(powerSpecMat, 2);
flux=specFlux(powerSpecMat, prm);

if prm.plot
	subplot(3,1,1);
	time=(1:length(wObj.signal))/wObj.fs;
	plot(time, wObj.signal);
	set(gca, 'xlim', [min(time), max(time)]);
%	xlabel('Time (sec)'); title(sprintf('Waveform of "%s"', strrep(waveFile, '_', '\_')));
	xlabel('Time (sec)'); title('Waveform');
	subplot(3,1,2);
	imagesc(powerSpecMat); axis xy
	title('Spectrogram');
	subplot(3,1,3);
	frameTime=frame2sampleIndex(1:frameNum-1, prm.frameSize, prm.overlap)/wObj.fs;
	plot(frameTime, flux, '.-');
	set(gca, 'xlim', [min(frameTime), max(frameTime)]);
	xlabel('Time (sec)'); ylabel('Flux');
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
flux=wave2specFlux(wObj, prm);
