function pitch = myFrame2pitch(frame, plotOpt, PP)

if nargin<1,    selfdemo;      return; end

% 求 frame 的平均值
average=mean(frame);
frame=frame-average;
maxShift=length(frame);
%method=1;
%plotOpt=1;
%acf=frame2acf(frame,maxShift,method);
%acf=myFrame2acf(frame);
%[pitch, index1, index2] = myAcf2pitch(acf, fs);
nsdf=myframe2nsdf(frame);
[pitch, index1, index2, clarity] = myNsdf2pitch(nsdf, PP.fs);
if frame2volume(frame)<PP.volTh
	%if plotOpt, fprintf('音量太小，音高設定為 0！\n'); end
	pitch=0;
end
if clarity<0.85
	%if plotOpt, fprintf('清澈度太小，音高設定為 0！\n'); end
	pitch=0;
end

if plotOpt

	subplot(3,1,1);
    plot(frame, '.-'); axis tight; title('Input frame');
  subplot(3,1,2);
    plot(nsdf, '.-'); title('nsdf'); axis tight;
  subplot(3,1,3);
    xVec=1:length(nsdf);
    localMaxIndex=[index1,index2];
    plot(xVec, nsdf, '.-', xVec(localMaxIndex), nsdf(localMaxIndex), 'ro');
    axis tight; title(sprintf('nsdf vector index1= %d  index2=%d   pitch= %4g   clarity=%g',index1,index2,pitch,clarity));
end

%=====selfdemo
function selfdemo
waveFile='soo.wav';
[y, fs, nbits]=wavread(waveFile);
y=y*2^nbits/2;
framedY=buffer(y, 256, 0);
frame=framedY(:, 220);
plotOpt=1;
method=1;
feval(mfilename, frame, fs, method, plotOpt);