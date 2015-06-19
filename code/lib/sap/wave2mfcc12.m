function mfcc=wave2mfcc12(y, fs, mfccOpt, plotOpt)
% wave2mfcc12: Conversion from wave to 12-dimensional mfcc matrix
%	Usage: mfcc=wave2mfcc12(y, mfccOpt, plotOpt)

%	Roger Jang, 20060421

if nargin<1, selfdemo; return; end
if nargin<2, fs=16000; end
if nargin<3, mfccOpt=mfccOptSet(fs); end
if nargin<4, plotOpt=0; end
PrmframeSize=mfccOpt.frameSize;
overlap=mfccOpt.overlap;

% Pre-emphasize
y2 = filter([1, -mfccOpt.preEmCoef], 1, y);
frameMat=buffer2(y2, frameSize, overlap);
frameNum=size(frameMat, 2);
frameTime=((0:frameNum-1)*(frameSize-overlap)+0.5*frameSize)/fs;
mfcc=[];
for i=1:frameNum
	frame=frameMat(:, i);
	mfcc=[mfcc, frame2mfcc(frame, fs, mfccOpt.tbfNum, mfccOpt.mfccNum)];
end

if plotOpt
	pcolor(mfcc); shading flat; axis xy
end

% ====== Selfdemo
function selfdemo
waveFile='主人下馬客在船.wav';
waveFile='Since his wife died, his life has had no meaning_70.wav';
[y, fs, nbits]=wavreadInt(waveFile);
mfccOpt=mfccOptSet(fs);
plotOpt=1;
mfcc=feval(mfilename, y, fs, mfccOpt, plotOpt);