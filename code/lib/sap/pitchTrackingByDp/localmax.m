function index = localmax(x, plotOpt)
% LOCALMAX Find local maxima of a vector
%	The output is a 0-1 index vector indicating where the local maxima are.

%	Roger Jang, 1999

if nargin==0, selfdemo; return; end
if nargin<2, plotOpt=0; end

b1 = x(2:end-1)>=x(1:end-2);
b2 = x(2:end-1)>x(3:end);
index = 0*x;
index(2:end-1) = b1 & b2;
if x(1)>x(2), index(1)=1; end
if x(end)>x(end-1), index(end)=1; end

if plotOpt,
	time = 1:length(x);
	plot(time, x, time(logical(index)), x(logical(index)), 'ro');
end

function selfdemo
t = 0:0.2:2*pi;
x = sin(t)+randn(size(t));
index = localmax(x, 1);