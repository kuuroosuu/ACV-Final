function feaMat = vec2frameFeaMat(vec, prm)
% vec2frameFeaMat: Vector to frame feature matrix
%	Usage: feaMat = vec2frameFeaMat(y, prm)
%
%	See selfdemo for a usage example.

%	Roger Jang, 20100130

if nargin<1, selfdemo; return; end
if nargin<2, prm=[]; end

prm=prmSet(prm);
frameSize=prm.frameSize;
overlap=prm.overlap;
frameFcn=prm.frameFcn;
frameFcnPrm=prm.frameFcnPrm;

vec=vec(:);
frameStep=frameSize-overlap;
frameCount=floor((length(vec)-overlap)/frameStep);

frameMat = zeros(frameSize, frameCount);
for i=1:frameCount
	startIndex=(i-1)*frameStep+1;
	frameMat(:, i)=vec(startIndex:(startIndex+frameSize-1));
end

if isempty(frameFcn)
	feaMat=frameMat;
	return;
end

feaMat=[];
for i=1:frameCount
	feaMat=[feaMat, feval(frameFcn, frameMat(:,i), frameFcnPrm)];
end

if prm.plot
	subplot(2,1,1);
	plot(vec);
	xlabel('Samples');
	subplot(2,1,2);
	imagesc(feaMat); axis xy
	title('Frame feature matrix');
	xlabel('Frame index');
end

% ====== Set default parameters
function prm=prmSet(prm)
if isempty(prm)
	clear prm
	prm.frameSize=256;
	prm.overlap=0;
	prm.frameFcn=[];
	prm.frameFcnPrm=[];
end
if ~isfield(prm, 'frameSize'), prm.frameSize=256; end
if ~isfield(prm, 'overlap'), prm.overlap=0; end
if ~isfield(prm, 'frameFcn'), prm.frameFcn=[]; end
if ~isfield(prm, 'frameFcnPrm'), prm.frameFcnPrm=[]; end

% ====== Self demo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
%waveFile='wubai_solicitude_orig.wav';
[y, fs, nbits]=wavReadInt(waveFile);
prm.frameSize=256;
prm.overlap=128;
prm.frameFcn='fftOneSide';
magSpecMat=vec2frameFeaMat(y, prm);
powerSpecMat=20*log10(magSpecMat);
freq=(0:prm.frameSize/2)*fs/prm.frameSize;
frameNum=size(magSpecMat, 2);
frameTime=frame2sampleIndex(1:frameNum, prm.frameSize, prm.overlap)/fs;