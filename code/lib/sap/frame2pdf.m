function [out, frameEstimated, excitation]=frame2pdf(frame, opt, showPlot)
% frame2acf: PDF (periodicity detection function) of a given frame (primarily for pitch tracking)
%
%	Usage:
%		out=frame2pdf(frame, opt, showPlot);
%			frame: Given frame
%			opt: Options for PDF computation
%				opt.pdf: PDF function to be used
%					'acf' for ACF
%					'amdf' for AMDF
%					'nsdf' for NSDF
%					'acfOverAmdf' for ACF divided by AMDF
%				opt.maxShift: no. of shift operations, which is equal to the length of the output vector
%				opt.method: 1 for using the whole frame for shifting
%					2 for using the whole frame for shifting, but normalize the sum by it's overlap area
%					3 for using frame(1:frameSize-maxShift) for shifting
%				opt.useSift: 1 for using SIFT
%				opt.SiftOrder: order of SIFT
%			showPlot: 0 for no plot, 1 for plotting the frame and ACF output
%			out: the returned PDF vector
%
%	Example:
%		waveFile='soo.wav';
%		wObj=waveFile2obj(waveFile);
%		frameSize=256;
%		frameMat=enframe(wObj.signal, frameSize);
%		frame=frameMat(:, 292);
%		opt=frame2pdf('defaultOpt');
%		subplot(4,1,1); plot(frame, '.-');
%		title('Input frame'); axis tight
%		subplot(4,1,2);
%		opt.method=1; opt.maxShift=length(frame);
%		out=frame2pdf(frame, opt);
%		plot(out, '.-'); title('method=1'); axis tight
%		subplot(4,1,3);
%		opt.method=2;
%		out=frame2pdf(frame, opt);
%		plot(out, '.-'); title('method=2'); axis tight
%		subplot(4,1,4);
%		opt.method=3;
%		opt.maxShift=length(frame)/2;
%		out=frame2pdf(frame, opt);
%		plot(out, '.-'); title('method=3'); axis tight
%
%	See also frame2amdf, frame2nsdf.

%	Roger Jang 20020404, 20041013, 20060313

if nargin<1, selfdemo; return; end
if nargin==1 && ischar(frame) && strcmpi(frame, 'defaultOpt')	% Set default options
	out.pdf='acf';
	out.maxShift=512;
	out.method=1;
	out.useSift=0;
	out.siftOrder=20;
	return
end
if nargin<2||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

frameEstimated=[];
excitation=[];
if opt.useSift
	[frameEstimated, excitation, coef]=sift(frame, opt.siftOrder);	% Simple inverse filtering tracking
	frame=excitation;
end
frameSize=length(frame);
maxShift=opt.maxShift;
out=zeros(maxShift, 1);

switch lower(opt.pdf)
	case 'acf'
	%	out=frame2acf(frame, maxShift, opt.method);
		out=frame2acfMex(frame, maxShift, opt.method);
	case 'amdf'
	%	out=frame2amdf(frame, maxShift, opt.method);
		out=frame2amdfMex(frame, maxShift, opt.method);
	case 'nsdf'
	%	out=frame2nsdf(frame, maxShift, opt.method);
		out=frame2nsdfMex(frame, maxShift, opt.method);
	case 'amdf4pt'
		opt.pdf='amdf';
		amdf=feval(mfilename, frame, opt);
		out=max(amdf)*(1-linspace(0,1,length(amdf))')-amdf;
	case 'acfoveramdf'
		opt.pdf='acf';
		acf =feval(mfilename, frame, opt);
		opt.pdf='amdf';
		amdf=feval(mfilename, frame, opt);
		out=0*acf;
		out(2:end)=acf(2:end)./amdf(2:end);
	otherwise
		error('Unknown PDF=%s!', opt.pdf);
end

if showPlot
	subplot(2,1,1);
	plot(frame, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('Input frame');
	subplot(2,1,2);
	plot(out, '.-');
	set(gca, 'xlim', [-inf inf]);
	title(sprintf('%s vector (opt.method = %d)', opt.pdf, opt.method));
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
