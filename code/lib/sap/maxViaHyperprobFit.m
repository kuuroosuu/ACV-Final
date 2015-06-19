function [optimPos, optimValue, coef]=optimViaHyperprobFit(vec, showPlot)
%optimViaHyperprobFit: Find the min/max of a equally spaced vector via hyperproble fitting
%	Usage:
%		[optimPos, optimValue, coef]=optimViaHyperprobFit(vec, showPlot)
%
%	Description:
%		[optmPos, optimValue, coef]=optimViaHyperprobFit(vec, showPlot)
%			vec: A given 3-element vector with the maximum/minimum at the central position
%
%	Example:
%		vec=randn(1,3);
%		[optimValue, optimIndex]=min(vec);
%		if optimIndex~=2, vec(optimIndex)=vec(2); vec(2)=optimValue; end
%		optimPos=optimViaHyperprobFit(vec, 1);

if nargin<1, selfdemo; return; end
if nargin<2, showPlot=1; end

c=vec(2);
b=(vec(3)-vec(1))/2;
a=(vec(3)+vec(1))/2-c;
optimPos=-b/(2*a);
optimValue=c-b*b/(4*a);
coef=[a b c];

if showPlot
	x=linspace(-1.5, 1.5);
	y=polyval(coef, x);
	plot(x, y, [-1 0 1], vec(:)', 'o');
	line(optimPos, optimValue, 'marker', 'o', 'color', 'r');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);