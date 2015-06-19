function feaVec = vibratoFeaExtract(pitchWindow, volumeWindow, showPlot)
% vibratoFeaExtract: Feature extraction for vibrato detection
%
%	Usage:
%		feaVec = vibratoFeaExtract(pitchWindow, volumeWindow)
%
%	Description:
%		feaVec = vibratoFeaExtract(pitchWindow, volumeWindow) returns the feature vector for vibrato detection
%
%	Example:
%		waveFile='vibrato02.wav';
%		wObj=waveFile2obj(waveFile);
%		pfType=1;		% 0 for AMDF, 1 for ACF
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
%		ptOpt.mainFun='maxPickingOverPf';	% This method can capture fine variations in pitch
%		ptOpt.useEpd=0;
%		ptOpt.useVolThreshold=0;		% volRatio is determine by epdByVol.
%		ptOpt.useClarityThreshold=0;
%		ptOpt.usePitchSmooth=0;
%		pitch=pitchTracking(wObj, ptOpt);
%		opt=wave2volume('defaultOpt');
%		opt.frameSize=ptOpt.frameSize; opt.overlap=ptOpt.overlap;
%		volume=wave2volume(wObj, opt);
%		vibratoFeaExtract(pitch(:), volume(:), 1);
%

if nargin<1, selfdemo; return; end
if nargin<3, showPlot=0; end

%% ====Pitch
y0 = pitchWindow(:);
time=(1:length(y0))';
poly=polyfit(time, y0, 2);
trend=polyval(poly, time);
y1=y0-trend;
m = fit(time, y1, 'sin1');
a = m.a1;	%振幅
b = m.b1;	%frequency
c = m.c1;	%shift
y2=m.a1*sin(m.b1*time+m.c1);
diffVec=y2-y1;
distance=norm(diffVec)/sqrt(length(diffVec));
%% ====Volume
x1 = 1:length(volumeWindow);
P = polyfit(x1(:),volumeWindow,2);%先做polyfit
for k = 1:length(volumeWindow)
    volumeWindow(k) = volumeWindow(k)-P(1)*k^2-P(2)*k-P(3);%減去二次趨勢
end
volume_window1 = volumeWindow-mean(volumeWindow);%做水平平移
m1 = fit(x1(:),volume_window1(:),'sin1');%做sin fit
v_a = m1.a1;
v_b = m1.b1;
v_c = m1.c1;
v_distance = 0;
for k = 1:length(volumeWindow) %caculate distance
	v_distance = v_distance + abs(volume_window1(k)-v_a*sin(v_b*k*1.0+v_c)-mean(volume_window1))^2;
end
v_distance = sqrt(v_distance);

feaVec=[a,b,distance,v_a,v_b,v_distance]';

if showPlot
	pitchLen=length(pitchWindow);
	time=1:pitchLen;
	time2=linspace(min(time), max(time));;
	pitch2=m.a1*sin(m.b1*time2+m.c1);
	subplot(211); plot(time, y0, time, y1, time2, pitch2);
	volLen=length(volumeWindow);
	time=1:volLen;
	vol2=m1.a1*sin(m1.b1*time+m1.c1)+polyval(P, time);
	subplot(212); plot(time, volumeWindow, time, vol2);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
