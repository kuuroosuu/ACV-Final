function [wObjOut2, wObjOut1]=audioCompress(wObj, opt, plotOpt)
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
	wObjOut2.threshold=-30;	% -10 dB
	wObjOut2.range=[-60, 0];
	wObjOut2.ratio=8;
	wObjOut2.smoothRange=6;	% 6dB starting from the knee
	wObjOut2.makeupGain=-wObjOut2.threshold*(1-1/wObjOut2.ratio);		% 8 dB by default
	wObjOut2.alpha1=0.001;	% 大==>爬升快 (For fs=44100)
	wObjOut2.alpha2=0.01;	% 大==>下降快 (For fs=44100)
%	wObjOut2.attackTime=1;	% ms
%	wObjOut2.releaseTime=10;	% ms
%	wObjOut2.releaseCheckTime=20;	% ms
	return
end
if nargin<2||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, plotOpt=0; end

if opt.threshold<opt.range(1) | opt.threshold>opt.range(2)
	error('The given threshold should be in the range of %s! (Current value=%g)', mat2str(opt.range), opt.threshold);
end

% Assume the starting point is [x0, y0] in the box, then 
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

x=wObj.signal;
o1=x;	% Output with time-decay effect
o2=x;	% Output without time-decay effect
dbo1=zeros(size(x));
dbo2=zeros(size(x));
gain1=zeros(size(x));
gain2=zeros(size(x));
signal1=zeros(size(x));
signal2=zeros(size(x));
p0=1;		% For computing decibels
for i=1:length(x)
	dbx=20*log10(abs(x(i))/p0);
	dby=staticTransferFunction(dbx, opt, x0, y0, x1, y1, a, b, c);
	gain1(i)=dby-dbx;
	if i>1
		if ~isnan(gain1(i))
			if gain1(i)>gain2(i-1)
				gain2(i)=opt.alpha1*gain1(i)+(1-opt.alpha1)*gain2(i-1);
			else
				gain2(i)=opt.alpha2*gain1(i)+(1-opt.alpha2)*gain2(i-1);
			end
		else
			gain2(i)=gain2(i-1);
		end
	end
	dbo1(i)=dbx+gain1(i);
	dbo2(i)=dbx+gain2(i);
	signal1(i)=p0*10^(dbo1(i)/20)*sign(x(i));
	signal2(i)=p0*10^(dbo2(i)/20)*sign(x(i));
	if signal2(i)>1, signal2(i)=1; end
	if signal2(i)<-1, signal2(i)=-1; end
end
wObjOut1=wObj; wObjOut1.signal=signal1; wObjOut1.file='';
wObjOut2=wObj; wObjOut2.signal=signal2; wObjOut2.file='';

if plotOpt==1
	time=1000*(1:length(wObj.signal))'/wObj.fs;
	linGain=wObjOut2.signal./wObj.signal;
	subplot(2,1,1);
	plot(time, wObj.signal, '-'); xlabel('Time (ms)'); grid on
	set(gca, 'ylim', [-1.1, 1.1]);
	wavePlayButton(wObj);
	subplot(2,1,2);
	plot(time, wObjOut2.signal, '-'); xlabel('Time (ms)'); grid on
	set(gca, 'ylim', [-1.1, 1.1]);
	wavePlayButton(wObjOut2);
	threshold=p0*10^(opt.threshold/20);
	subplot(2,1,1);
	axisLimit=axis;
	line([axisLimit(1:2), axisLimit(2), nan, axisLimit(1), axisLimit(1:2)], [-threshold*[1 1], -threshold, nan, threshold, threshold*[1 1]], 'color', 'r');
	line(time, wObjOut1.signal./wObj.signal, 'color', 'g');
	line(time, wObjOut2.signal./wObj.signal, 'color', 'm');
	legend('Waveform', 'Threshold', 'Ideal linear gain', 'Linear gain', 'location', 'northOutside', 'orientation', 'horizontal');
	subplot(2,1,2);
	axisLimit=axis;
	line([axisLimit(1:2), axisLimit(2), nan, axisLimit(1), axisLimit(1:2)], [-threshold*[1 1], -threshold, nan, threshold, threshold*[1 1]], 'color', 'r');
	line(time, wObjOut1.signal./wObj.signal, 'color', 'g');
	line(time, wObjOut2.signal./wObj.signal, 'color', 'm');
	legend('Waveform', 'Threshold', 'Ideal linear gain', 'Linear gain', 'location', 'northOutside', 'orientation', 'horizontal');
end

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
		dby(i)=staticTransferFunction(dbx(i), opt, x0, y0, x1, y1, a, b, c);
	end
	line(dbx, dby, 'marker', '.', 'linestyle', 'none');
	subplot(122);
	x=linspace(eps, 1, 1001);
	for i=1:length(x)
		temp=staticTransferFunction(20*log10(x(i)/p0), opt, x0, y0, x1, y1, a, b, c);
		y(i)=p0*10^(temp/20);
	end
%	keyboard
	plot(x, y, '.-'); grid on; axis image; axis([0 1 0 1]);
	xlabel('Input (linear)'); ylabel('Output (linear)');
	title('I/O function in linear scale');
end


% ====== Single sample conversion
function output=staticTransferFunction(input, opt, x0, y0, x1, y1, a, b, c)
if isinf(input)		% To avoid the output being nan
	output=-inf;
	return
end
if input<x0
	dby=input;
elseif input<x1
	dby=a*input^2+b*input+c;
else
	dby=(input-x1)/opt.ratio+y1;
end
output=dby+opt.makeupGain;		% Use makeupGain

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
