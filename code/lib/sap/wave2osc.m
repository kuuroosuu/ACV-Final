function [oscObj, rawOsc, smoothedOsc, localMeanOsc, magSpec]=wave2osc(wObj, oscOpt, showPlot)
%wave2osc: Wave to onset strength curve (OSC) conversion
%
%	Usage:
%		[oscObj, rawOsc, smoothedOsc, localMeanOsc]=wave2osc(wObj, oscOpt, showPlot)
%			oscObj.signal: signal for the onset strength curve
%			oscObj.time: time for the onset strength curve
%
%	Example:
%		waveFile='song01s5.wav';
%		wObj=waveFile2obj(waveFile);
%		oscOpt=wave2osc('defaultOpt');
%		showPlot=1;
%		oscObj=wave2osc(wObj, oscOpt, showPlot);
%		wavePlayButton(wObj);	% audioPlayWithBar is called within this function!

%	Roger Jang, 20120330

if nargin<1, selfdemo; return; end
% ====== Set the default options
if nargin==1 && ischar(wObj) && strcmpi(wObj, 'defaultOpt')
	oscObj.fsTarget=8192;
	oscObj.frameSize=256;
	oscObj.frameStep=32;
	oscObj.numMelBin=40;	% Number of mel bins
	oscObj.divisor1=50;
	oscObj.divisor2=300;
	return
end
if nargin<2|isempty(oscOpt), oscOpt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

if isstr(wObj), wObj=waveFile2obj(wObj); end
wObj.signal=mean(wObj.signal, 2);	% Stereo ==> Mono
if (wObj.fs~=oscOpt.fsTarget)		% Resample to have the right sampling rate
	gg=gcd(wObj.fs, oscOpt.fsTarget);
	wObj.signal=resample(wObj.signal, oscOpt.fsTarget/gg, wObj.fs/gg);	% up-sampling¡GoscOpt.fsTarget/gg ­¿ down-sampling¡GwObj.fs/gg ­¿
	wObj.fs=oscOpt.fsTarget;
end

frameSize=oscOpt.frameSize;
frameStep=oscOpt.frameStep;
overlap=frameSize-frameStep;
numMelBin=oscOpt.numMelBin;		% Number of mel channels
wave=wObj.signal;
fs=wObj.fs;
frameRate=fs/frameStep;			% No of frames per second

% ====== Construct mel-scale spectrogram
[spec, binFreq, oscTime]=spectrogram(wave, frameSize, overlap, frameSize, fs);	% Compute spectrogram
magSpec=abs(spec);	% Magnitude spectrum
mlmx = fft2melmx(frameSize, fs, numMelBin);	% A weight matrix to combine FFT bins into mel bins (This can be computed in advance.)
magSpec = mlmx(:,1:(frameSize/2+1))*magSpec;	% Mel-scale amplitude spectrum

% ====== Delete some frequency bands if necessary
magSpec2=magSpec;
if isfield(oscOpt, 'deleteIndex')
	magSpec(oscOpt.deleteIndex, :)=nan;	% For plotting
	magSpec2(oscOpt.deleteIndex, :)=[];	% For computing OSC
end

% ====== Compute onset strength curve
rawOsc=mean(max(0, diff(magSpec2, 1, 2)));	% Raw onset strength curve after half-wave rectifier
oscTime(end)=[];				% Due to the above diff operation
%[meanSubtractedOsc, smoothedOsc, localMeanOsc] = gaussian_smooth_localmean(rawOsc, frameRate, frameRate/50, true, frameRate/300);
coef1=gaussianFilterCoef(frameRate, frameRate/oscOpt.divisor1);
smoothedOsc=gaussianFilter(rawOsc, frameRate, coef1);		% Smoothed OSC
coef2=gaussianFilterCoef(frameRate, frameRate/oscOpt.divisor2);
localMeanOsc=gaussianFilter(smoothedOsc, frameRate, coef2);	% Trend of OSC
meanSubtractedOsc=smoothedOsc-localMeanOsc;
meanSubtractedOsc(meanSubtractedOsc<0)=0;			% Trend-subtracted OSC
oscObj=struct('signal', meanSubtractedOsc, 'time', oscTime);

if showPlot
	clf;
	minTime=0; maxTime=(length(wObj.signal)-1)/wObj.fs;
	subplot(411);
	imagesc(oscTime, 1:size(magSpec,1), 20*log10(magSpec)); axis xy;
%	colormap(gray);
%	colorbar('location','northoutside');
	xlabel('Time (sec)'); ylabel('Freq. bin index'); title('Spectrogram');
	subplot(412);
	plot(oscTime, rawOsc, oscTime, smoothedOsc); set(gca, 'xlim', [minTime, maxTime]);
	xlabel('Time (sec)'); ylabel('Amplitude'); title('OSC (original and smoothed)');
	subplot(413);
	plot(oscTime, smoothedOsc); set(gca, 'xlim', [minTime, maxTime]);
	line(oscTime, localMeanOsc, 'color', 'r');
	xlabel('Time (sec)'); ylabel('Amplitude'); title('Smoothed OSC and its trend');
	subplot(414);
	plot(oscTime, meanSubtractedOsc); set(gca, 'xlim', [minTime, maxTime]);
	xlabel('Time (sec)'); ylabel('Amplitude'); title('Trend-subtracted OSC');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);