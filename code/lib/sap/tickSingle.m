function tick=tickSingle(opt, showPlot)
% tickSingle: Generate a single tick sound
%
%	Usage:
%		tick=tickSingle(opt, showPlot)
%
%	Example:
%		opt=tickSingle('defaultOpt');
%		tick=tickSingle(opt, 1);
%		sound(tick, opt.fs);

%	Roger Jang, 20130623

if nargin<1, selfdemo; return; end
% ====== Set the default options
if nargin==1 && ischar(opt) && strcmpi(opt, 'defaultOpt')
	tick.freq=600;
	tick.fs=16000;
	tick.duration=0.1;
	return
end
if nargin<2, showPlot=0; end

time=(0:opt.fs*opt.duration-1)/opt.fs;
y=sin(2*pi*opt.freq*time);
opt=noteEnvelope('defaultOpt');
opt.peakTime=0.003;
opt.width=0.005;
opt.max=1;
envelope=noteEnvelope(time, opt);
tick=y.*envelope;

if showPlot
	subplot(211);
	plot(time, y);
	subplot(212);
	plot(time, tick);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
