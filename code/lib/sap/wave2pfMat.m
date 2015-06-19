function pfMat=wave2pfMat(y, frameSize, overlap, pfLen, pfType, pfMethod, plotOpt)
% wave2pfMat: Conversion from wave to PF (pitch function, either AMDF or ACF) matrix
%	Usage: pfMat=wave2pfMat(y, frameSize, overlap, pfLen, pfType, pfMethod, plotOpt)
%		pfLen: PF length, or equivalently, max time delay for computing PF
%		pfType=0 for AMDF, 1 for ACF
%		pfMethod=1 (moving base = whole frame)
%			 2 (moving base = whole frame, but normalized by the overlap area)
%			 3 (moving base = frame(1:frameSize-pfLen))
%			 4 (moving base = frame(1:pfLen), assuming pfLen<=frameSize/2)

%	Roger Jang, 20050227, 20100606

if nargin<1, selfdemo; return; end
if nargin<2, frameSize=256; end
if nargin<3, overlap=0; end
if nargin<4, pfLen=frameSize/2; end
if nargin<5, pfType=0; end		% 0 for AMDF, 1 for ACF
if nargin<6, pfMethod=3; end
if nargin<7, plotOpt=0; end

pfMat=wave2pfMatMex(y, frameSize, overlap, pfLen, pfType, pfMethod);

if plotOpt
	pcolor(pfMat); shading flat; axis xy
end

% ====== Selfdemo
function selfdemo
waveFile='Since his wife died, his life has had no meaning_70.wav';
[y, fs, nbits]=wavReadInt(waveFile);
ptOpt=ptOptSet(fs, nbits);
pfLen=ptOpt.frameSize/2;
pfType=0;
if pfType==0, pfMethod=2; end
if pfType==1, pfMethod=1; end
plotOpt=1;
pfMat=feval(mfilename, y, ptOpt.frameSize, ptOpt.overlap, pfLen, pfType, pfMethod, plotOpt);
title(strPurify4label(['pfMat of "', waveFile, '"']));
