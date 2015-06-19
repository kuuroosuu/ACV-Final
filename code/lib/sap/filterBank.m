function [filterCoef, centerFreq]=filterBank(binNum, fs, filterOrder, plotOpt)
%filterBank: Coefficients for a filter bank
%	Usage: filterCoef=filterBank(binNum, fs, filterOrder, plotOpt)

%	Roger Jang, 20101213

if nargin<1, selfdemo; return; end
if nargin<2, fs=44100; end
if nargin<3, filterOrder=5; end
if nargin<4, plotOpt=0; end

maxFreq=fs/2;
halfSpan=maxFreq/binNum/2;
boundaryFreq=linspace(0, maxFreq, binNum+1);
centerFreq=(boundaryFreq(1:end-1)+boundaryFreq(2:end))/2;
for i=1:binNum
	passBand=(centerFreq(i)+halfSpan*[-1, 1])/maxFreq;
%	fprintf('%d/%d: passBand=%s\n', i, binNum, mat2str(passBand));
	if i==1
		[b, a]=butter(filterOrder, passBand(2), 'low');
	elseif i==binNum
		[b, a]=butter(filterOrder, passBand(1), 'high');
	else
		[b, a]=butter(filterOrder, passBand);
	end
	[h(:,i), w(:,i)]=freqz(b, a);
	filterCoef{i,1}=b;
	filterCoef{i,2}=a;
end

if plotOpt
	plot(w/pi*maxFreq, abs(h));
	set(gca', 'xlim', [0, maxFreq]);
	xlabel('Freq (Hz)');
	title('Filter bank'); grid on
end

% ====== Self demo
function selfdemo
binNum=20;
fs=44100;
filterOrder=5;
filterCoef=filterBank(binNum, fs, filterOrder, 1);
