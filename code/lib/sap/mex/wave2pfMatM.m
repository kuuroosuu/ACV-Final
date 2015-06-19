function pfMat=wave2pfMat(y, frameSize, overlap, maxShift, pfType, pfMethod, plotOpt)
% wave2pfMat: Conversion from wave to amdf matrix
%	Usage: pfMat=wave2pfMat(y, frameSize, overlap, maxShift, pfType, pfMethod, plotOpt)
%		maxShift: time delay for computing PF
%		pfType=0 for AMDF, 1 for ACF
%		pfMethod=1 (moving base = whole frame)
%			 2 (moving base = whole frame, but normalized by the overlap area)
%			 3 (moving base = frame(1:frameSize-maxShift))
%			 4 (moving base = frame(1:maxShift), assuming maxShift<=frameSize/2)

%	Roger Jang, 20050227

if nargin<1, selfdemo; return; end
if nargin<2, frameSize=256; end
if nargin<3, overlap=0; end
if nargin<4, maxShift=frameSize/2; end
if nargin<5, pfType=0; end		% ACF
if nargin<6, pfMethod=3; end
if nargin<7, plotOpt=0; end

frameMat=buffer2(y, frameSize, overlap);
frameNum=size(frameMat, 2);
pfMat=zeros(maxShift, frameNum);
switch(pfType)
	case 0	% AMDF
		for i=1:frameNum
			frame=frameMat(:, i);
			frame2=frameFlip(frame);
		%	frame2=localAverage(frame2);			% C program does not do this.
			pfMat(:,i)=frame2amdf(frame2, maxShift, pfMethod);
		end
	case 1	% ACF
		for i=1:frameNum
			frame=frameMat(:, i);
%			if i==1, asciiWrite(frame, 'frame0.txt'); end
			frame2=frameFlip(frame);
		%	frame2=frame2-fix(mean(frame2));		% 配合Ｃ，所以加上 fix (目前無法改成mean，因為C的STL的問題尚未解決...)
			frame2=frame2-mean(frame2);
%			if i==1, asciiWrite(frame2, 'frame1.txt'); end
		%	frame2=localAverage(frame2);			% C program does not do this.			
			pfMat(:,i)=frame2acf(frame2, maxShift, pfMethod);
		end
end

if plotOpt
	pcolor(pfMat); shading flat; axis xy
end

% ====== Selfdemo
function selfdemo
waveFile='主人下馬客在船.wav';
waveFile='Since his wife died, his life has had no meaning_70.wav';
[y, fs, nbits]=wavreadInt(waveFile);
PP=ptOptSet(fs, nbits);
maxShift=PP.frameSize/2;
pfType=0;
pfMethod=3;
plotOpt=1;
pfMat=feval(mfilename, y, PP.frameSize, PP.overlap, maxShift, pfType, pfMethod, plotOpt);
title(['pfMat of "', waveFile, '"']);