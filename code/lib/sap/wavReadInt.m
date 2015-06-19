function [y, fs, nbits, opts, cueLabel]=wavReadInt(waveFile, range, plotOpt)
% waveReadInt: Same as wavread, but return integer value of y
%	Usage: [y, fs, nbits, opts, cueLabel]=wavReadInt(waveFile, range, plotOpt)
%		This first 4 output arguments are the same as those of builtin function wavread.m, except that
%		the values of y are integers. More specifically:
%			If nbits=8, then -128<=y<=127
%			If nbits=16, then -32768<=y<=32767
%		The fifth output cueLabel is the cue labels that were added by CoolEdit.
%	
%	For example:
%		waveFile='tapping.wav';
%		[y, fs, nbits, opts, cueLabel] = wavReadInt(waveFile);
%		time=((1:length(y))/fs);
%		plot(time, y); set(gca, 'xlim', [-inf inf]);
%		axisLimit=axis;
%		for i=1:length(cueLabel)
%			line(cueLabel(i)*[1, 1]/fs, axisLimit(3:4), 'color', 'r'); 
%		end

%	Roger Jang, 20050228, 20090607

if nargin==0, selfdemo; return; end
if nargin==1, range=inf; plotOpt=0; end
if nargin==2
	if length(range)==1, plotOpt=range; range=inf; end	% waveReadInt(waveFile, plotOpt)
	if length(range)==2, plotOpt=0; end			% waveReadInt(waveFile, range)
end

if isempty(range) | isinf(range)
	[y, fs, nbits, opts]=wavread(waveFile);
else
	[y, fs, nbits, opts]=wavread(waveFile, range);
end
y=y*2^nbits/2;

if nargout>4
	fid=fopen(waveFile, 'rb');
	waveContents = fread(fid);
	fclose(fid);
	cueLabel=binWave2cueLabel(waveContents);
end

if plotOpt
	channelNum=size(y, 2);
	time=((1:length(y))/fs);
	if channelNum==1
		plot(time, y); axis([-inf inf 2^nbits/2*[-1 1]]);
		axisLimit=axis;
		if nargout>4
			for i=1:length(cueLabel)
				line(cueLabel(i)*[1, 1]/fs, axisLimit(3:4), 'color', 'r'); 
			end
		end
	end
	if channelNum==2
		subplot(2,1,1); plot(time, y(:,1)); axis([-inf inf 2^nbits/2*[-1 1]]);
		axisLimit=axis;
		if nargout>4
			for i=1:length(cueLabel)
				line(cueLabel(i)*[1, 1]/fs, axisLimit(3:4), 'color', 'r'); 
			end
		end
		subplot(2,1,2); plot(time, y(:,2)); axis([-inf inf 2^nbits/2*[-1 1]]);
		axisLimit=axis;
		if nargout>4
			for i=1:length(cueLabel)
				line(cueLabel(i)*[1, 1]/fs, axisLimit(3:4), 'color', 'r'); 
			end
		end
	end
end

% ====== Self demo
function selfdemo
waveFile='tapping.wav';
plotOpt=1;
[y, fs, nbits, opts, cueLabel] = feval(mfilename, waveFile, plotOpt);