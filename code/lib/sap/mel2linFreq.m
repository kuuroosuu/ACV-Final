function linFreq=mel2linFreq(melFreq)
%mel2linFreq: Mel frequency to linear frequency conversoin

%	Roger Jang, 20020502

if nargin==0; selfdemo; return; end
linFreq=700*(exp(melFreq/1125)-1);			% melFreq=1125*ln(1+linFreq/700)

% ====== Self demo
function selfdemo
melFreq=0:3000;
linFreq=feval(mfilename, melFreq);
plot(melFreq, linFreq);
xlabel('Mel-frequency');
ylabel('Frequency');
title('Mel-frequency to frequency curve');
axis equal tight