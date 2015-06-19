function specDiffOut = specSub(y, opt, plotOpt)
% specDiffOut = specSub(y, opt)  
%
if nargin<1, selfdemo; return; end
if nargin<2
	opt.frameSize=1024;
	opt.overlap=1024-128;
	opt.diffType=1;			% 1: keep positive, 2: abs
end
if nargin<3, plotOpt=0; end

frameMat1=buffer2(y(:,1), opt.frameSize, opt.overlap);
frameMat2=buffer2(y(:,2), opt.frameSize, opt.overlap);
frameNum=size(frameMat1, 2);
frameMat1=frameMat1.*(hamming(opt.frameSize)*ones(1,frameNum));
frameMat2=frameMat2.*(hamming(opt.frameSize)*ones(1,frameNum));
% ====== FFT
spec1=fft(frameMat1); spec1=spec1(1:opt.frameSize/2+1, :); mag1=spec1.*conj(spec1);
spec2=fft(frameMat2); spec2=spec2(1:opt.frameSize/2+1, :); mag2=spec2.*conj(spec2);
specDiff=mag1-mag2;

if opt.diffType==1
	specDiff=(specDiff+abs(specDiff))/2;	% Keep only the positive part
else
	specDiff=abs(specDiff);
end

%specDiffOut=0;
specDiffOut=istft(specDiff.^(1/2).*exp(1i*angle(spec1)), opt.frameSize, opt.frameSize,  opt.frameSize-opt.overlap)';

if plotOpt
	subplot(3,1,1); imagesc(db(mag1)); axis xy
	subplot(3,1,2); imagesc(db(mag2)); axis xy
	subplot(3,1,3); imagesc(db(specDiff)); axis xy
	xlabel('Frame index');
	ylabel('Frequency bin index');
end

% ====== Self demo
function selfdemo
demoFile='demo_stereo.wav';
[y, fs, nbits]=wavread(demoFile);
opt.frameSize=1024;
opt.overlap=1024-128;
opt.diffType=1;		% 1: positive  other:abs
plotOpt=1;
specDiffOut=specSub(y, opt, plotOpt);
