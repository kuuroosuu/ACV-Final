function flux=specNovelty(powerSpecMat, noveltyPrm, plotOpt)
% specFlux: Compute the flux of power spectrum
%	Usage: flux=specFlux(powerSpecMat, freq, freqRange)
%		powerSpecMat: matrix of power spectrum (in linear scale), where
%			each column corresponds to the power spectrum of a frame
%		noveltyPrm: paramters for spectral flux
%			noveltyPrm.useSpecNormalize: use spectrum normalization or not
%			noveltyPrm.useHalfRectifier: use half rectifier or not
%				(This usually set to 1 for onset detection.)
%		flux: output flux
%
%	For usage example, see selfdemo() of this function.

%	Roger Jang, 20070506, 20080401, 20100131

if nargin<1, selfdemo; return; end
if nargin<2 || isempty(noveltyPrm)
	noveltyPrm.useSpecNormalize=0;
	noveltyPrm.useHalfRectifier=1;
	noveltyPrm.filterOrder=7;
end
if nargin<3, plotOpt=0; end

flux=specFlux(powerSpecMat, noveltyPrm, plotOpt);
trend=meanFilter(flux, noveltyPrm.filterOrder);
novelty=flux-trend;
novelty(novelty<0)=0;

if plotOpt
	line(1:length(trend), trend, 'color', 'r');
	line(1:length(novelty), novelty, 'color', 'm', 'linewidth', 2);
end
legend('Flux', 'Trend', 'Novelty');

% ====== Self demo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
waveFile='song01.wav';
[y, fs, nbits]=wavread(waveFile);
y=y(4*fs:10*fs);
prm.frameSize=32/1000*fs;		% frameSize = 32 ms
prm.overlap=prm.frameSize-10/1000*fs;	% frameStep = 10 ms ===> 100 frames/sec
prm.frameFcn='fftOneSide';
prm.frameFcnPrm.fs=fs;
prm.frameFcnPrm.nfft=4*prm.frameSize;	% Zero-padding to this length
magSpecMat=vec2frameFeaMat(y, prm);
powerSpecMat=20*log10(magSpecMat);
plotOpt=1;
novelty=specNovelty(powerSpecMat, [], plotOpt);
