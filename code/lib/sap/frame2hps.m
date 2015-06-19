function [hps, freq, spec0, spec1, spec2]=frame2hps(frame, fs, zeroPaddedFactor, showPlot)
% frame2hps: Conversion from a frame to harmonic product spectrum, for pitch tracking
%
%	Usage: hps=frame2hps(frame, fs, zeroPaddedFactor, showPlot)
%
%	Example:
%		waveFile='soo.wav';
%		[y, fs, nbits]=wavread(waveFile);
%		startIndex=15000;
%		frameSize=round(32*fs/1000);
%		zeroPaddedFactor=15;
%		frame=y(startIndex:startIndex+frameSize-1);
%		[hps, freq]=frame2hps(frame, fs, zeroPaddedFactor, 1);

if nargin<1; selfdemo; return; end
if nargin<3, zeroPaddedFactor=15; end
if nargin<4; showPlot=0; end

frame2=frame.*hamming(length(frame));	% Hamming windowing
frameSize=length(frame);
frame2=[frame2; zeros(zeroPaddedFactor*frameSize, 1)];		% zero padding
[magSpec, phase, freq, powerSpec]=fftOneSide(frame2, fs, 1);
spec0=powerSpec;

% === zero mean
polyOrder=10;
m=length(spec0);
xPos=(1:m)'/m;
specTrend=polyval(polyfit(xPos, spec0, polyOrder), xPos);
spec1=spec0-specTrend;
% === Weighting
h=0.999;
spec2=spec1.*(h.^(0:length(spec1)-1)');

N=8;
component=nan*ones(length(spec2), N);
for i=1:N
	toBeAdded=spec2(1:i:end);
	component(1:length(toBeAdded), i)=toBeAdded;
end
hps=nanmean(component, 2);
%hps=nansum(component, 2);

if showPlot
	plotNum=5;
	subplot(plotNum,1,1);
	plot(1:length(frame), frame, '.-'); axis tight;
	title('Frame'); xlabel('Samples');

	subplot(plotNum,1,2);
	plot(freq, spec0, freq, specTrend);
	set(gca, 'xlim', [-inf inf]);
	title('Power spectrum and its trend');
	
	subplot(plotNum,1,3);
	plot(freq, spec1, freq, spec2);
	set(gca, 'xlim', [-inf inf]);
	title('Trend-subtracted power spectrum and its tapering version');
	
	subplot(plotNum,1,4);
	temp=component;
%	temp(temp==0)=nan;
	plot(freq, temp);
	set(gca, 'xlim', [-inf inf]);
	title('Down-sampled versions of power spectrum');

	subplot(plotNum,1,5);
	plot(freq, hps);
	set(gca, 'xlim', [-inf inf]);
	axisLimit=axis;
	title('Harmonic product spectrum');
	xlabel('Freq (Hz)');
	line(40*[1 1], axisLimit(3:4), 'color', 'r');
	line(1000*[1 1], axisLimit(3:4), 'color', 'r');
	hps2=hps; hps2(freq<40|freq>1000)=-inf;
	[maxValue, maxIndex]=max(hps2);
	pitch=freq2pitch(freq(maxIndex));
	line(freq(maxIndex), hps2(maxIndex), 'marker', 'o', 'color', 'r');
	fprintf('Pitch = %f semitone\n', pitch);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
