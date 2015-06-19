function ashod = frame2ashod(frameMat, order);
% frame2volume: Frame (or frame matrix) to absolute sum of high-order difference (ASHOD)
%
%	Usage:
%		volume = frame2ashod(frameMat)
%		volume = frame2ashod(frameMat, order)
%
%	Description:
%		volume = frame2ashod(frameMat, order)
%			frameMat: a column vector of a frame, or a matrix where each column is a frame
%			order: order of difference
%
%	Example:
%		waveFile='清華大學資訊系.wav';
%		waveFile='star_noisy.wav';
%		[y, fs, nbits]=wavReadInt(waveFile);
%		frameSize=256; overlap=0;
%		frameMat=buffer2(y, frameSize, overlap);
%		frameNum=size(frameMat, 2);
%		ashod1=frame2ashod(frameMat, 1);
%		ashod2=frame2ashod(frameMat, 2);
%		ashod3=frame2ashod(frameMat, 3);
%		ashod4=frame2ashod(frameMat, 4);
%		time=(1:length(y))/fs;
%		frameTime=((0:frameNum-1)'*(frameSize-overlap)+0.5*frameSize)/fs;
%		subplot(2,1,1); plot(time, y); axis tight; ylabel('Amplitude'); title(strPurify4label(waveFile));
%		subplot(2,1,2); plot(frameTime, [ashod1; ashod2; ashod3; ashod4]', '.-'); axis tight; ylabel('ASHOD');
%		legend('order=1', 'order=2', 'order=3', 'order=4');

%	Roger Jang, 20070417

if nargin<1, selfdemo; return; end
if nargin<2, order=4; end

diffMat=diff(frameMat, order);
ashod=sum(abs(diffMat));

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
