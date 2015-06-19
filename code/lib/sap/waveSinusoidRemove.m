function waveSinusoidRemove(inWavFile, outWavFile, freq, plotOpt)

if nargin<1, selfdemo; return; end
if nargin<3, freq=50; end
if nargin<4, plotOpt=0; end

[y, fs, nbits]=wavread(inWavFile);
%y=y(end-5000:end);
time=(1:length(y))'/fs;

% Linear search over freq+[-0.5, 0.5]
freqs=linspace(freq-0.5, freq+0.5, 101);
for i=1:length(freqs)
	f=freqs(i);
	A=[ones(length(y),1), cos(2*pi*f*time), sin(2*pi*f*time)];
	x=A\y;
	y2=A*x;
	diff=y-y2;
	error(i)=sum(diff.^2)/length(diff);
end
[minError, minIndex]=min(error);
theFreq=freqs(minIndex);
A=[ones(length(y),1), cos(2*pi*theFreq*time), sin(2*pi*theFreq*time)];
x=A\y;
y2=A*x;

fprintf('Writing output to %s...\n', outWavFile);
wavwrite(y2, fs, nbits, outWavFile);

if plotOpt
	fprintf('Best freq = %f\n', theFreq);
	subplot(3,1,1); plot(freqs, error, 'b.-'); line(freqs(minIndex), error(minIndex), 'marker', 'o', 'color', 'r'); xlabel('Freq'); ylabel('Error');
	subplot(3,1,2); plot(time, y, time, y2, 'r'); title('y'); z=axis;
	subplot(3,1,3); plot(time, y-y2); title('y-y3'); axis(z);
end

% ====== Self demo
function selfdemo
inFile='speechWith50HzNoise02.wav';
outFile=[tempname, '.wav'];
freq=50;
plotOpt=1;
feval(mfilename, inFile, outFile, freq, plotOpt);