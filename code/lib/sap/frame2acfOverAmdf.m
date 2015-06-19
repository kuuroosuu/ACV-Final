function [acfOverAmdf, acf, amdf]=frame2acfOverAmdf(frame, maxShift, method, showPlot)
% frame2acfOverAmdf: ACF/AMDF of a given frame (primarily for pitch tracking)
%
%	Usage:
%		out=frame2acfOverAMDF(frame, maxShift, method, showPlot);
%			maxShift: no. of shift operations, which is equal to the length of the output vector
%			method: 1 for using the whole frame for shifting
%				2 for using the whole frame for shifting, but normalize the sum by it's overlap area
%				3 for using frame(1:frameSize-maxShift) for shifting
%			showPlot: 0 for no plot, 1 for plotting the frame and ACF output
%			out: the returned ACF/AMDF vector
%
%	Example:
%		waveFile='soo.wav';
%		wObj=waveFile2obj(waveFile);
%		frameSize=256;
%		frameMat=enframe(wObj.signal, frameSize);
%		frame=frameMat(:, 292);
%		subplot(4,1,1); plot(frame, '.-');
%		title('Input frame'); axis tight
%		subplot(4,1,2);
%		method=1; out=frame2acfOverAmdf(frame, frameSize, method);
%		plot(out, '.-'); title('method=1'); axis tight
%		subplot(4,1,3);
%		method=2; out=frame2acfOverAmdf(frame, frameSize, method);
%		plot(out, '.-'); title('method=2'); axis tight
%		subplot(4,1,4);
%		method=3; out=frame2acfOverAmdf(frame, round(frameSize/2), method);
%		plot(out, '.-'); title('method=3'); axis tight
%
%	See also frame2acf, frame2amdf, frame2nsdf.

%	Roger Jang 20020404, 20041013, 20060313

if nargin<1, selfdemo; return; end
if nargin<2, maxShift=floor(length(frame)/2); end
if nargin<3, method=3; end
if nargin<4, showPlot=0; end

acf =frame2acf (frame, maxShift, method);
amdf=frame2amdf(frame, maxShift, method);
acfOverAmdf=0*acf;
acfOverAmdf(2:end)=acf(2:end)./amdf(2:end);

if showPlot
	plotNum=4;
	subplot(plotNum,1,1);
	set(plot(frame, '.-'), 'tag', 'frame');
	axis tight;
	title('Frame');
	subplot(plotNum,1,2);
	plot(acf, '.-');
	axis tight;
	title(sprintf('ACF vector (method = %d)', method));
	subplot(plotNum,1,3);
	plot(amdf, '.-');
	axis tight;
	title(sprintf('AMDF vector (method = %d)', method));
	subplot(plotNum,1,4);
	plot(acfOverAmdf, '.-');
	axis tight
	title(sprintf('ACF/AMDF vector (method = %d)', method));
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);