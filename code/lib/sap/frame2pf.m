function pf=frame2pf(frame, pfType, maxShift, method, plotOpt)
% frame2pf: COmpute PF (pitch function, including ACF, AMDF, NSDF, etc) from a given frame
%
%	Usage:
%		out=frame2pf(frame, pfType, maxShift, method, plotOpt);
%			frame: input frame
%			pfType: 0 for AMDF, 1 for ACF, 2 for NSDF
%			maxShift: no. of shift operations, which is equal to the length of the output vector
%			method: 1 for using the whole frame for shifting
%				2 for using the whole frame for shifting, but normalize the sum by it's overlap area
%				3 for using frame(1:frameSize-maxShift) for shifting
%				4 for using frame(1:maxShift) for shifting
%			plotOpt: 0 for no plot, 1 for plotting the frame and ACF output
%			pf: the returned PF vector
%
%	Example:
%		waveFile='greenOil.wav';
%		[y, fs, nbits]=wavReadInt(waveFile);
%		frameMat=buffer2(y, 256, 0);
%		frame=frameMat(:, 250);
%		subplot(4,1,1); plot(frame, '.-');
%		title('Input frame'); axis tight
%		subplot(4,1,2); pfType=0; out=frame2pf(frame, pfType);
%		plot(out, '.-'); title('pfType=0'); axis tight
%		subplot(4,1,3); pfType=1; out=frame2pf(frame, pfType);
%		plot(out, '.-'); title('pfType=1'); axis tight
%		subplot(4,1,4); pfType=2; out=frame2pf(frame, pfType);
%		plot(out, '.-'); title('pfType=2'); axis tight
%
%	See also frame2acf, frame2amdf, frame2nsdf.

%	Roger Jang 20020404, 20041013, 20060313, 20100617

if nargin<1, selfdemo; return; end
if nargin<2, pfType=1; end
if nargin<3, maxShift=length(frame); end
if nargin<4, method=1; end
if nargin<5, plotOpt=0; end

switch pfType
	case 0
		pf=frame2amdf(frame, maxShift, method, plotOpt);
	case 1
		pf=frame2acf(frame, maxShift, method, plotOpt);
	case 2
		pf=frame2nsdf(frame, maxShift, method, plotOpt);
	otherwise
		error('Unknown method!');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);