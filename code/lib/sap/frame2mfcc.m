function mfcc=frame2mfcc(frame, fs, filterNum, mfccNum, plotOpt)
% frame2mfcc: Frame to MFCC conversion.
%	Usage: mfcc=frame2mfcc(frame, fs, filterNum, mfccNum, plotOpt)
%
%	For example:
%		waveFile='what_movies_have_you_seen_recently.wav';
%		[y, fs, nbits]=wavReadInt(waveFile);
%		startIndex=12000;
%		frameSize=512;
%		frame=y(startIndex:startIndex+frameSize-1);
%		frame2mfcc(frame, fs, 20, 12, 1);

%	Roger Jang 20060417

if nargin<1, selfdemo; return; end
if nargin<2, fs=16000; end
if nargin<3, filterNum=20; end
if nargin<4, mfccNum=12; end
if nargin<5, plotOpt=0; end

frameSize=length(frame);
% ====== Preemphasis should be done at wave level
%a=0.95;
%frame2 = filter([1, -a], 1, frame);
frame2=frame;
% ====== Hamming windowing
frame3=frame2.*hamming(frameSize);
% ====== FFT
[fftMag, fftPhase, fftFreq, fftPowerDb]=fftOneSide(frame3, fs);
% ====== Triangular band-pass filter bank
triFilterBankPrm=getTriFilterBankPrm(fs, filterNum);	% Get parameters for triangular band-pass filter bank
% Triangular bandpass filter.
for i=1:filterNum
	tbfCoef(i)=dot(fftPowerDb, trimf(fftFreq, triFilterBankPrm(:,i)));
end
% ====== DCT
mfcc=zeros(mfccNum, 1);
for i=1:mfccNum
	coef = cos((pi/filterNum)*i*((1:filterNum)-0.5))';
	mfcc(i) = sum(coef.*tbfCoef');
end
% ====== Log energy
%logEnergy=10*log10(sum(frame.*frame));
%mfcc=[logEnergy; mfcc];

if plotOpt
	subplot(2,1,1);
	plot(frame, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('Input frame');
	subplot(2,1,2);
	plot(mfcc, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('MFCC vector');
end

% ====== trimf.m (from fuzzy toolbox)
function y = trimf(x, prm)
a = prm(1); b = prm(2); c = prm(3);
y = zeros(size(x));
% Left and right shoulders (y = 0)
index = find(x <= a | c <= x);
y(index) = zeros(size(index));
% Left slope
if (a ~= b)
    index = find(a < x & x < b);
    y(index) = (x(index)-a)/(b-a);
end
% right slope
if (b ~= c)
    index = find(b < x & x < c);
    y(index) = (c-x(index))/(c-b);
end
% Center (y = 1)
index = find(x == b);
y(index) = ones(size(index));

% ====== Self demo
function selfdemo
waveFile='what_movies_have_you_seen_recently.wav';
[y, fs, nbits]=wavReadInt(waveFile);
startIndex=12000;
frameSize=512;
frame=y(startIndex:startIndex+frameSize-1);
feval(mfilename, frame, fs, 20, 12, 1);