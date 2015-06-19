function waveHighPassFilter(inWavFile, outWavFile, cutOffFreq, filterOrder, plotOpt)

if nargin<1, selfdemo; return; end
if nargin<3, cutOffFreq=100; end
if nargin<4, filterOrder=2; end
if nargin<5, plotOpt=0; end

[y, fs, nbits]=wavread(inWavFile);
%y=y(end-2000:end);
time=(1:length(y))'/fs;

[b, a]=butter(filterOrder, cutOffFreq/(fs/2), 'high');
y2=filter(b, a, y);

minIndex=1;
% Find the best shift
%for shift=1:200
%	difference=y2(shift:end)-y(1:end-shift+1);
%	error(shift)=sum(difference.^2)/length(difference);
%end
%[minError, minIndex]=min(error);
%y2=y2(minIndex:end);
fprintf('Writing output to %s...\n', outWavFile);
wavwrite(y2, fs, nbits, outWavFile);

if plotOpt
	fprintf('minIndex=%d\n', minIndex);
%	plot(error);
%	line(minIndex, error(minIndex), 'marker', 'o', 'color', 'r');
	% ====== Plot the result
	subplot(3,1,1);
	plot(time, y, 'b', time, y2, 'r'); title('y and y2'); z=axis;
	subplot(3,1,2);
	y3=[y2(minIndex:end); nan*ones(minIndex-1, 1)];
	plot(time, y, 'b', time, y3, 'r'); title('y and y3'); z=axis;
	subplot(3,1,3);
	plot(time, y-y3, 'b'); title('y-y3');
end

% ====== Self demo
function selfdemo
inWavFile='speechWith50HzNoise02.wav';
outWavFile=[tempname, '.wav'];
feval(mfilename, inWavFile, outWavFile, 100, 5, 1);