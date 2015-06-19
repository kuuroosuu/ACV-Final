function y=tapDelayLine(x, prm, plotOpt)
% tapDelayLine: Tap delay line

if nargin<1, selfdemo; return; end
if nargin<2|isempty(prm), prm=tapDelayLinePrmSet; end
if nargin<3, plotOpt=0; end

totalDelay=sum(prm.delay);
delayLine=zeros(1, totalDelay);	% memory allocation for length 10
nextPos=0;
for i=1:length(x);
	nextPos=nextPos+1;
	if nextPos>totalDelay, nextPos=1; end
	y(i)=0;
	delayIndex=nextPos;
	for j=1:length(prm.delay)
		delayIndex=delayIndex-prm.delay(j);
		if delayIndex<1, delayIndex=delayIndex+totalDelay; end
		y(i)=y(i)+prm.gain(j)*delayLine(delayIndex);
	end	
	delayLine(nextPos)=x(i);
end

if plotOpt
	subplot(2,1,1); stem(x, 'filled'); title('Input');
	subplot(2,1,2); stem(y, 'filled'); title('Output');
end

% ====== Self demo
function selfdemo
x=zeros(100,1); x(1)=1;		% unit impulse signal of length 100
prm.delay=[10 24 37];
prm.gain=[0.8 0.6 0.3];
plotOpt=1;
y=tapDelayLine(x, prm, plotOpt);
