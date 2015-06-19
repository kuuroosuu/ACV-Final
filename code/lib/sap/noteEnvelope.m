function out = noteEnvelope(timeVec, opt, showPlot)
% noteEnvelope: Envelope of a note
%
%	Usage:
%		out = noteEnvelope(timeVec, opt, showPlot)
%	
%	Description:
%		out = noteEnvelope(timeVec, opt) returns a note envelope with parameters a and b.
%
%	Example:
%		% === Plot all three envelope functions
%		opt=noteEnvelope('defaultOpt');
%		time=linspace(0, 0.3, 1001);
%		showOpt=1;
%		subplot(311); envelope1=noteEnvelope(time, opt, 1); title(opt.type);
%		opt.type='exponential';
%		subplot(312); envelope1=noteEnvelope(time, opt, 1); title(opt.type);
%		opt.type='exponential2';
%		subplot(313); envelope1=noteEnvelope(time, opt, 1); title(opt.type);
%		% === Plot of the family of 'simple' envelope function
%		opt.type='simple';
%		opt.width=0.02;
%		time=linspace(0, 0.3, 1001);
%		envelope1=noteEnvelope(time, opt);
%		opt.width=0.04;
%		envelope2=noteEnvelope(time, opt);
%		opt.width=0.06;
%		envelope3=noteEnvelope(time, opt);
%		figure; plot(time, envelope1, time, envelope2, time, envelope3);

%	Roger Jang, 20130414

if nargin<1, selfdemo; return; end
if nargin==1 && ischar(timeVec) && strcmpi(timeVec, 'defaultOpt')	% Set default options
	out.type='simple';
	out.max=0.9;			% Max. amplitude of the sound wave
	out.peakTime=0.015;		% Parameters for 'simple' envelope
	out.width=0.05;			% Parameters for 'simple' envelope
	out.k=10;		% For 'exponential' envelope
	out.period=0.05;
	out.exponent=20;
	return
end
if nargin<2||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

switch lower(opt.type)
	case 'simple'
		% The function = c*t/(t^2+a*t+b^2), with the peak at [b, c/(a+2*b)]
		% Let z=a+4*b, then 50% height occurs at (z-sqrt(z*z-4*b^2))/2 and (z+sqrt(z*z-4*b^2))/2.
		b=opt.peakTime;
		a=sqrt(opt.width^2+4*b^2)-4*b;
		c=opt.max*(a+2*b);
	%	fprintf('a=%g, b=%g, c=%g\n', a, b, c);
		out = c*timeVec./(timeVec.^2+a*timeVec+b^2);
	case 'exponential'
		out = 1-exp((timeVec-max(timeVec))*opt.k);
		out=out/(max(out))*opt.max;
	case 'exponential2';
		out = sin(2*pi/opt.period*timeVec).*exp(-opt.exponent*timeVec);
end

if showPlot
	plot(timeVec, out);
%	line(sqrt(b), opt.max, 'color', 'r', 'marker', 'o');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
