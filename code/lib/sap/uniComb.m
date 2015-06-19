function y=uniComb(x, prm, plotOpt)
% uniComb: Universal comb filter

if nargin<1, selfdemo; return; end
if nargin<2|isempty(prm), prm=uniCombPrmSet; end
if nargin<3, plotOpt=0; end

delayLine=zeros(1, prm.delay);	% memory allocation for length 10
delayIndex=1;
for i=1:length(x);
	xh=x(i)+prm.fb*delayLine(delayIndex);
	y(i)=prm.ff*delayLine(delayIndex)+prm.bl*xh;
	delayLine(delayIndex)=xh;
	delayIndex=delayIndex+1;
	if delayIndex>prm.delay, delayIndex=1; end
end;

if plotOpt
	subplot(2,1,1); stem(x, 'filled'); title('Input');
	subplot(2,1,2); stem(y, 'filled'); title('Output');
end

% ====== Self demo
function selfdemo
x=zeros(100,1); x(1)=1;		% unit impulse signal of length 100
prm.bl=0.5;
prm.fb=-0.5;
prm.ff=1;
prm.delay=10;
plotOpt=1;
y=uniComb(x, prm, plotOpt);
