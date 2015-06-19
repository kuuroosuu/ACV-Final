function y=lpComb(x, prm, plotOpt)
% lpComb: Low-pass comb filter

if nargin<1, selfdemo; return; end
if nargin<2|isempty(prm), prm=lpCombPrmSet; end
if nargin<3, plotOpt=0; end

delayedX=zeros(1, prm.delay);
delayedZ=zeros(1, prm.delay);
delayIndex=1;
for i=1:length(x)
	index1=delayIndex-prm.delay; if (index1<1) index1=index1+prm.delay; end
	index2=delayIndex-1; if (index2<1) index2=index2+prm.delay; end
%	fprintf('i=%d, index1=%d, index2=%d\n', i, index1, index2);
	z=delayedX(index1)+prm.g*delayedZ(index1)+prm.a*delayedZ(index2);
	y(i)=delayedX(index1)+prm.g*delayedZ(index1);
	delayedX(delayIndex)=x(i);
	delayedZ(delayIndex)=z;
	delayIndex=delayIndex+1;
	if (delayIndex>prm.delay) delayIndex=1; end
end;

if plotOpt
	subplot(2,1,1); stem(x, 'filled'); title('Input');
	subplot(2,1,2); stem(y, 'filled'); title('Output');
end

% ====== Self demo
function selfdemo
x=zeros(100,1); x(1)=1;		% unit impulse signal of length 100
prm.delay=10;
prm.a=0.2;
prm.g=0.5;
plotOpt=1;
y=lpComb(x, prm, plotOpt);
