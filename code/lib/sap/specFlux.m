function flux=specFlux(powerSpecMat, sfPrm)
% specFlux: Compute the flux of power spectrum
%	Usage: flux=specFlux(powerSpecMat, freq, freqRange)
%		powerSpecMat: matrix of power spectrum (in linear scale), where
%			each column corresponds to the power spectrum of a frame
%		sfPrm: paramters for spectral flux
%			sfPrm.useSpecNormalize: use spectrum normalization or not
%			sfPrm.useHalfRectifier: use half rectifier or not
%				(This usually set to 1 for onset detection.)
%		flux: output flux
%
%	For usage example, see selfdemo() of this function.

%	Roger Jang, 20070506, 20080401, 20100131

if nargin<1, selfdemo; return; end
if nargin<2, sfPrm=[]; end
sfPrm=prmSet(sfPrm);

% If there is -inf in powerSpecMat
[p, q]=find(isinf(powerSpecMat));
if ~isempty(p)
	temp=powerSpecMat(:);
	temp(find(isinf(temp)))=[];
	minValue=min(temp);
	for i=1:length(p)
		powerSpecMat(p(i), q(i))=minValue;
	end
end

if sfPrm.useSpecNormalize
	for i=1:size(powerSpecMat, 2)
		minValue=min(powerSpecMat(:,i));
		maxValue=max(powerSpecMat(:,i));
		powerSpecMat(:,i)=(powerSpecMat(:,i)-minValue)/(maxValue-minValue);
	end
end

diffMat=diff(powerSpecMat, 1, 2);
if sfPrm.useHalfRectifier
	diffMat=(diffMat+abs(diffMat))/2;
end

flux=sum(abs(diffMat));

if sfPrm.plot
	subplot(2,1,1);
	imagesc(powerSpecMat); axis xy
	title('Spectrogram');
	subplot(2,1,2);
	plot(flux, '.-');
	xlabel('Frame index'); ylabel('Flux');
end

% ====== Set default parameters
function prm=prmSet(prm)
% ====== Default prm
def.useSpecNormalize=0;
def.useHalfRectifier=1;
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
[y, fs, nbits]=wavread(waveFile);
prm.frameSize=32/1000*fs;		% frameSize = 32 ms
prm.overlap=prm.frameSize-10/1000*fs;	% frameStep = 10 ms ===> 100 frames/sec
prm.frameFcn='fftOneSide';
prm.frameFcnPrm.fs=fs;
prm.frameFcnPrm.nfft=4*prm.frameSize;	% Zero-padding to this length
magSpecMat=vec2frameFeaMat(y, prm);
powerSpecMat=20*log10(magSpecMat);
sfPrm.plot=1;
flux=specFlux(powerSpecMat, sfPrm);
