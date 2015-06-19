function out=frame2nsdf(frame, maxShift, method, plotOpt)
% frame2nsdf: NSDF (normalized sum of difference function) of a given frame (primarily for pitch tracking)
%
%	Usage:
%		out=frame2nsdf(frame, maxShift, method, plotOpt);
%			maxShift: no. of shift operations, which is equal to the length of the output vector
%			method: 1 for using the whole frame for shifting
%				2 for using the whole frame for shifting, but normalize the sum by it's overlap area
%				3 for using frame(1:frameSize-maxShift) for shifting
%			plotOpt: 0 for no plot, 1 for plotting the frame and ACF output
%			out: the returned ACF vector
%
%	Example:
%		waveFile='greenOil.wav';
%		wObj=waveFile2obj(waveFile);
%		frameSize=256;
%		frameMat=enframe(wObj.signal, frameSize);
%		frame=frameMat(:, 250);
%		subplot(4,1,1);
%		plot(frame, '.-');
%		title('Input frame'); axis tight
%		subplot(4,1,2);
%		method=1; out=frame2nsdf(frame, 256, method);
%		plot(out, '.-'); title('method=1'); axis tight
%		subplot(4,1,3);
%		method=2; out=frame2nsdf(frame, 256, method);
%		plot(out, '.-'); title('method=2'); axis tight
%		subplot(4,1,4);
%		method=3; out=frame2nsdf(frame, 128, method);
%		plot(out, '.-'); title('method=3'); axis tight
%
%	Reference: "A SMARTER WAY TO FIND PITCH" by Philip McLeod, Geoff Wyvill
%
%	See also frame2acf, frame2amdf.

%	Roger Jang 20070209

if nargin<1, selfdemo; return; end
if nargin<2, maxShift=length(frame); end
if nargin<3, method=1; end
if nargin<4, plotOpt=0; end

frameSize=length(frame);
out=zeros(maxShift, 1);
switch method
	case 1		% moving base = whole frame
		for i=1:maxShift
			out(i) = dot(frame(1:frameSize-i+1), frame(i:frameSize));
			deminator=sum(frame(1:frameSize-i+1).^2)+sum(frame(i:frameSize).^2);
			if deminator>0, out(i) = 2*out(i)/deminator; end
		end
	case 2		% moving base = whole frame, but normalized by the overlap area
		for i=1:maxShift
			out(i) = dot(frame(1:frameSize-i+1), frame(i:frameSize));
			deminator=sum(frame(1:frameSize-i+1).^2)+sum(frame(i:frameSize).^2);
			if deminator>0, out(i) = 2*out(i)/deminator; end
		end
	case 3		% moving base = frame(1:frameSize-maxShift)
		for i=1:maxShift
			out(i) = dot(frame(1:frameSize-maxShift), frame(i:frameSize-maxShift+i-1));
			deminator=sum(frame(1:frameSize-maxShift).^2)+sum(frame(i:frameSize-maxShift+i-1).^2);
			if deminator>0, out(i) = 2*out(i)/deminator; end
		end
	otherwise
		error('Unknown method!');
end

if plotOpt
	subplot(2,1,1);
	plot(frame, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('Input frame');
	subplot(2,1,2);
	plot(out, '.-');
	set(gca, 'xlim', [-inf inf]);
	title(sprintf('NSDF vector (method = %d)', method));
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
