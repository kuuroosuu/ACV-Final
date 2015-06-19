function [optimPos, optimValue, coef]=optimViaParabolaFit(vec, showPlot)
%optimViaParabolaFit
%optimViaParabolaFit: Find the min/max position based on a given 3-element vector
%
%	Usage:
%		[optimPos, optimValue, coef]=optimViaParabolaFit(vec, showPlot)
%
%	Example:
%		% === Find the max position
%		subplot(1,2,1);
%		vec=rand(1,3);
%		[maxValue, maxIndex]=max(vec);
%		if maxIndex~=2, vec(maxIndex)=vec(2); vec(2)=maxValue; end
%		optimViaParabolaFit(vec, 1);
%		% === Find the min position
%		subplot(1,2,2);
%		vec=rand(1,3);
%		[maxValue, maxIndex]=min(vec);
%		if maxIndex~=2, vec(maxIndex)=vec(2); vec(2)=maxValue; end
%		optimViaParabolaFit(vec, 1);

%	Roger Jang, 20121107

if nargin<1, selfdemo; return; end
if nargin<2, showPlot=0; end

if vec(1)==vec(2) && vec(2)==vec(3)
	optimPos=0;
	optimValue=vec(1);
	coef=[0 0 vec(1)];
else
	c=vec(2);
	b=(vec(3)-vec(1))/2;
	a=(vec(3)+vec(1))/2-c;
	optimPos=-b/(2*a);
	optimValue=c-b*b/(4*a);
	coef=[a b c];	
end

if showPlot
	x=linspace(-1.5, 1.5);
	y=polyval(coef, x);
	plot(x, y, [-1 0 1], vec(:)', 'o');
	line(optimPos, optimValue, 'marker', '^', 'color', 'r');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);