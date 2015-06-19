function melFreq=lin2melFreq(linFreq)
%lin2melFreq: Linear frequency to mel frequency conversion

%	Roger Jang, 20020502

if nargin==0; selfdemo; return; end
melFreq=1125*log(1+linFreq/700);

% ====== Self demo
function selfdemo
linFreq=0:8000;
melFreq=feval(mfilename, linFreq);
plot(linFreq, melFreq);
xlabel('Frequency');
ylabel('Mel-frequency');
title('Frequency to mel-frequency curve');
axis equal tight