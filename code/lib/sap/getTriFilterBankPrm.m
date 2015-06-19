function freq=getTriFilterBankPrm(fs, filterNum, plotOpt)
% getTriFilterBankPrm: Get the parameters of the triangular band-pass filter bank used in computing MFCC
%	Usage: freq=getTriFilterBankPrm(fs, filterNum, plotOpt)

%	Roger Jang, 20060417

if nargin<1; selfdemo; return; end
if nargin<2, filterNum=20; end
if nargin<3, plotOpt=0; end

fLow=0;
fHigh=fs/2;
% Compute the frequencies of the triangular band-pass filters
for i=1:filterNum+2
	f(i)=mel2linFreq(lin2melFreq(fLow)+(i-1)*(lin2melFreq(fHigh)-lin2melFreq(fLow))/(filterNum+1));
end
freq=[];
for i=1:filterNum
	freq=[[f(i); f(i+1); f(i+2)], freq];
end

if plotOpt==1
	% Plot the triangular band-pass filters
	filter=[];
	for i=1:filterNum
		filter=[[0; 1; 0], filter];
	end
	subplot(2,1,1);
	plot(freq, filter);
	xlabel('Frequency (Hz)');
	title('Triangular filter bank');
	% Plot the normalized triangular band-pass filters
	filter=[];
	for i=1:filterNum
		filter=[[0; 2/(f(i+2)-f(i)); 0], filter];
	end
	subplot(2,1,2);
	plot(freq, filter);
	xlabel('Frequency (Hz)');
	title('Triangular filter bank (normalized)');
end

% ====== Self demo
function selfdemo()
fs=16000;
filterNum=20;
plotOpt=1;
feval(mfilename, fs, filterNum, plotOpt);