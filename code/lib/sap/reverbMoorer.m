function y3=reverbMoorer(x, prm, plotOpt);

if nargin<1, selfdemo; return; end
if nargin<2 | isempty(prm), prm=reverbMoorerPrmSet; end
if nargin<3, plotOpt=0; end

% ====== Step 1: Tap delay line (pre-delay)
tdlPrm.delay=round(prm(1:18)); prm(1:18)=[];
tdlPrm.gain=prm(1:18); prm(1:18)=[];
t0=clock;
y1=tapDelayLineMex(x, tdlPrm);
%fprintf('tap delay line ===> %.2f sec\n', etime(clock, t0));

% ====== Step 2: 6 parallel uniComb/lpComb filter
lpCombPrmDelay=round(prm(1:6)); prm(1:6)=[];
lpCombPrmA=prm(1:6); prm(1:6)=[];
lpCombPrmG=prm(1:6); prm(1:6)=[];
output=zeros(length(x), length(lpCombPrmDelay));
t0=clock;
for i=1:length(lpCombPrmDelay)
	lpCombPrm=lpCombPrmSet;
	lpCombPrm.delay=lpCombPrmDelay(i);
	lpCombPrm.a=lpCombPrmA(i);
	lpCombPrm.g=lpCombPrmG(i);
	output(:,i)=lpCombMex(y1, lpCombPrm);
end
y2=sum(output, 2);
y2=y2-mean(y2);
y2=y2/max(abs(y2));
%fprintf('6 comb filters ===> %.2f sec\n', etime(clock, t0));

% ====== Step 3: Final all-pass filter
apfPrm.delay=round(prm(1)); prm(1)=[];
apfPrm.bl=prm(1); prm(1)=[];
apfPrm.fb=-apfPrm.bl;
apfPrm.ff=1;
t0=clock;
y3=uniCombMex(y2, apfPrm);
%y3=uniCombMex(y3, apfPrm);
%fprintf('All-pass filter ===> %.2f sec\n', etime(clock, t0));

if plotOpt
	subplot(4,1,1); stem(x, 'filled');
	subplot(4,1,2); stem(y1, 'filled');
	subplot(4,1,3); stem(y2, 'filled');
	subplot(4,1,4); stem(y3, 'filled');
end

% ====== Self demo
function selfdemo
fs=44100; x=zeros(2*fs, 1); x(1)=1;	% This is for obtaining the impulse response
plotOpt=1;
y=reverbMoorer(x, [], plotOpt);