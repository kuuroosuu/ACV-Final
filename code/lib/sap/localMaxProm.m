function controlPoints=localMaxProm(vec, opt, plotOpt)
% localMaxProm: Get prominent local maxima
%
%	Usage:
%		controlPoints=localMaxProm(vec, opt, plotOpt)
%
%	Example:
%		time=1:500;
%		vec=abs(sin(2*pi*0.002*time))+randn(1,500)/10;
%		opt.thresholdRatio=0.05;
%		opt.minSeparation=10;
%		controlPoints=localMaxProm(vec, opt, 1);

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

vecTh=(max(vec)-min(vec))*opt.thresholdRatio+min(vec);
onset1=vec>vecTh & localMax(vec);

onset2=onset1;
halfWinWidth=opt.minSeparation;
index=find(onset2);
for i=index
	startIndex=max(i-round(halfWinWidth), 1);
	endIndex=min(i+round(halfWinWidth), length(vec));
	[junk, maxIndex]=max(vec(startIndex:endIndex));
	if maxIndex+startIndex-1~=i
		onset2(i)=0;
	end
end
controlPoints=onset2;

if plotOpt
	plot(1:length(vec), vec, '.-'); set(gca, 'xlim', [-inf inf]);
	xlabel('Sample index'); ylabel('Vec');
	line([1, length(vec)], vecTh*[1 1], 'color', 'r');
	line(find(onset1), vec(onset1), 'marker', '.', 'color', 'g', 'linestyle', 'none');
	line(find(onset2), vec(onset2), 'marker', '^', 'color', 'k', 'linestyle', 'none');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
