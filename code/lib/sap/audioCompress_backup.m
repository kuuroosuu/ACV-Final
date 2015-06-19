function wObj2=audioCompress(wObj, opt, plotOpt)
% audioCompress: Audio sample compressor
%
%	Usage:
%		wObj2=audioCompress(wObj)
%		wObj2=audioCompress(wObj, opt)
%		wObj2=audioCompress(wObj, opt, plotOpt)
%		audioCompress([], opt, 'transferFunctionPlot')
%
%	Description:
%		wObj2=audioCompress(wObj, opt) returns the output after audio compressor
%			wObj: input audio object
%			opt: options for the compressor ([] for using the defaults)
%				The default options can be obtained via opt=audioCompress('defaultOpt')
%			wObj2: output audio object
%		audioCompress([], opt, 'transferFunctionPlot') gives the plot of the static transfer function of the compressor.
%
%	Example:
%		% === Try audio compressor on a speech utterance
%		wObj=waveFile2obj('what_movies_have_you_seen_recently.wav');
%		wObj.signal=wObj.signal(1:1.0*wObj.fs);
%		%wObj.fs=16000; time=(1:0.5*wObj.fs)/wObj.fs; wObj.signal=sin(2*pi*200*time);
%		%index=1:0.1*wObj.fs;	wObj.signal(index)=wObj.signal(index)*0.29;
%		%index=0.3*wObj.fs:0.5*wObj.fs;	wObj.signal(index)=wObj.signal(index)*0.29;
%		opt=audioCompress('defaultOpt');
%		wObj2=audioCompress(wObj, opt, 1);
%		% === Plot the static I/O characteristics
%		figure; audioCompress([], opt, 'transferFunctionPlot');
%
%	Roger Jang, 20120704

if nargin<1, selfdemo; return; end
if ischar(wObj) && strcmpi(wObj, 'defaultOpt')	% Set the default options
	wObj2.threshold=-10;	% -10 dB
	wObj2.range=[-40, 0];
	wObj2.ratio=5;
	wObj2.smoothRange=6;	% 6dB starting from the knee
	wObj2.makeupGain=8;		% 8 dB by default
	wObj2.attackTime=10;	% ms
	wObj2.releaseTime=40;	% ms
	wObj2.releaseCheckTime=20;	% ms
	return
end
if nargin<2||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, plotOpt=0; end

if opt.threshold<opt.range(1) | opt.threshold>opt.range(2)
	error('The given threshold should be in the range of %s! (Current value=%g)', mat2str(opt.range), opt.threshold);
end

% Assume the starting point is [x0, x0] in the box, then 
% The curve is y = ax^2+(1-2*a*x0)*x+a*x0^2
knee=opt.range(2)+opt.threshold;
x0=knee-opt.smoothRange;
y0=x0;
x1=knee+opt.smoothRange;
y1=knee+opt.smoothRange/opt.ratio;
if y1>x1
	y1=x1;	% This happens when opt.ratio<1
	x1=knee+opt.smoothRange*opt.ratio;
end
a=(y1-x1)/(x1^2-2*x0*x1+x0^2);
b=1-2*a*x0;
c=a*x0^2;

if isempty(wObj), wObj.signal=[]; wObj.fs=0; end

signal=wObj.signal;

signal2=signal;
signChange=0;
attackStartIndex=0;
attackStartAcDegree=0;
releaseStartIndex=-inf;
prevOverIndex=0;
acDegree=zeros(1, length(signal));	% Overall AC degree, used in attack and release phases
status=zeros(1, length(signal));	% Status of the compression: 0 for no compression, 1 for attack, 2 for full compression, 3 for release
prevStatus=0;
p0=1;		% For computing decibels
for i=1:length(signal)
	x=signal(i);
	if x<0, x=-x; signChange=1; end
	dbx=20*log10(x/p0);
	if dbx>opt.threshold, prevOverIndex=i; end	% For release check time
	status(i)=prevStatus;	% Default current state is the previous state
	switch(prevStatus)
		case 0		% No compression
			if dbx>=opt.threshold
				status(i)=1;
				attackStartIndex=i;
			end
			acDegree(i)=0;
		case 1		% Attack phase (where compression degree is increasing)
			timeSinceAttack=(i-attackStartIndex)/wObj.fs*1000;
			if timeSinceAttack>=opt.attackTime;
				status(i)=2;
			end
			acDegree(i)=attackStartAcDegree+(1-attackStartAcDegree)*timeSinceAttack/opt.attackTime;	% Overall AC degree due to attack
		case 2		% Full compression
			if (i-prevOverIndex)/wObj.fs*1000>opt.releaseCheckTime
				status(i)=3;
				releaseStartIndex=i;
			end
			acDegree(i)=1;	% Ratio
		case 3		% Release phase (where compression degree is decreasing)
			timeSinceRelease=(i-releaseStartIndex)/wObj.fs*1000;			
			acDegree(i)=max(0, 1-timeSinceRelease/opt.releaseTime);	% Ratio
			if timeSinceRelease>opt.releaseTime
				status(i)=0;
			end
			if dbx>opt.threshold
				status(i)=1;
				attackStartIndex=i;
				attackStartAcDegree=acDegree(i);
			end
		otherwise
			error('Unknown status!');
	end
	dby=singleSampleConversion(dbx, opt, x0, y0, x1, y1, a, b, c, acDegree(i));
	signal2(i)=p0*10^(dby/20);
%	if isinf(dbx), keyboard; end
	if signChange, signal2(i)=-signal2(i); end
	signChange=0;
	prevStatus=status(i);
end

wObj2=wObj; wObj2.signal=signal2; wObj2.file='';

if strcmp(plotOpt, 'transferFunctionPlot')
	subplot(121);
	plot(opt.range, opt.range); axis image; grid on
	title('I/O function in dB scale');
	line(opt.range, opt.threshold*[1 1], 'color', 'r');
	line([knee, opt.range(2)], [opt.range(2)+opt.threshold, opt.range(2)+opt.threshold+(opt.range(2)-knee)/opt.ratio], 'color', 'k');
	bb=[x0, y0, opt.smoothRange*2, opt.smoothRange*2];
	boxOverlay(bb, 'g', 1);
	line(x0, y0, 'marker', 'o', 'color', 'm');
	line(x1, y1, 'marker', 'o', 'color', 'm');
	x=linspace(x0, x1);
	y=a*x.^2+b*x+c;
	line(x, y, 'color', 'k');
	xlabel('Input (dB)');
	ylabel('Output (dB)');
%	x=linspace(eps, 1); dbx=20*log10(x/p0);
	dbx=linspace(opt.range(1), opt.range(2));
	dby=dbx;
	for i=1:length(dbx)
		dby(i)=singleSampleConversion(dbx(i), opt, x0, y0, x1, y1, a, b, c);
	end
	line(dbx, dby, 'marker', '.', 'linestyle', 'none');
	subplot(122);
	x=linspace(eps, 1, 1001);
	for i=1:length(x)
		temp=singleSampleConversion(20*log10(x(i)/p0), opt, x0, y0, x1, y1, a, b, c);
		y(i)=p0*10^(temp/20);
	end
%	keyboard
	plot(x, y, '.-'); grid on; axis image; axis([0 1 0 1]);
	xlabel('Input (linear)'); ylabel('Output (linear)');
	title('I/O function in linear scale');
	return
end

if plotOpt
	time=1000*(1:length(wObj.signal))'/wObj.fs;
	linGain=wObj2.signal./wObj.signal;
	subplot(2,1,1);
	plot(time, wObj.signal, '-'); xlabel('Time (ms)');
	set(gca, 'ylim', [-1.1, 1.1]);
	wavePlayButton(wObj);
	subplot(2,1,2);
	plot(time, wObj2.signal, '-'); xlabel('Time (ms)');
	set(gca, 'ylim', [-1.1, 1.1]);
	wavePlayButton(wObj2);
	threshold=p0*10^(opt.threshold/20);
	subplot(2,1,1);
	axisLimit=axis;
	line([axisLimit(1:2), axisLimit(2), nan, axisLimit(1), axisLimit(1:2)], [-threshold*[1 1], -threshold, nan, threshold, threshold*[1 1]], 'color', 'r');
	line(time, acDegree, 'color', 'm');
	line(time, status/4, 'color', 'g', 'linewidth', 1);
	line(time, acDegree, 'color', 'r');
	legend('Waveform', 'Threshold', 'Degree of compression', 'Status of compressor', 'Linear gain', 'location', 'northOutside', 'orientation', 'horizontal');
	subplot(2,1,2);
	axisLimit=axis;
	line([axisLimit(1:2), axisLimit(2), nan, axisLimit(1), axisLimit(1:2)], [-threshold*[1 1], -threshold, nan, threshold, threshold*[1 1]], 'color', 'r');
	line(time, acDegree, 'color', 'm');
	line(time, status/4, 'color', 'g', 'linewidth', 1);
	line(time, acDegree, 'color', 'r');
	legend('Waveform', 'Threshold', 'Degree of compression', 'Status of compressor', 'Linear gain', 'location', 'northOutside', 'orientation', 'horizontal');
end

% ====== Single sample conversion
function output=singleSampleConversion(input, opt, x0, y0, x1, y1, a, b, c, acDegree)
if isinf(input) && input<0	% To avoid the output being nan
	output=-inf;
	return
end
if nargin<10, acDegree=1; end
if input<x0
	dby=input;
elseif input<x1
	dby=a*input^2+b*input+c;
else
	dby=(input-x1)/opt.ratio+y1;
end
output=dby+opt.makeupGain;		% Use makeupGain
output=(output-input)*acDegree+input;	% Use AC degree during attack/release

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
