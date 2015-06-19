function waveCurve4cooledit(wObj, curve, frameRate, outWaveFile)
% waveCurve4cooledit: Create wave and its line (e.g., pitch or novelty curve) for visualization in cooledit.
%
%	Usage:
%		waveCurve4cooledit(wObj, curve, frameRate)
%		waveCurve4cooledit(wObj, curve, frameRate, outWaveFile)
%
%	Description:
%		waveCurve4cooledit(wObj, curve, frameRate) create a temp wave file with the left channel as the wave and the right channel as the curve (e.g., pitch or novelty) that are easily visualized in cooledit.
%			wObj: A give wave object
%			curve: The corresponding curve for wave, such as pitch or  novelty curve
%			frameRate: Frame rate of the given curve
%		waveCurve4cooledit(wObj, curve, frameRate, waveFile) saves the result to waveFile.
%
%	Example:
%		waveFile='vibrato.wav';
%		wObj=waveFile2obj(waveFile);
%		pfType=1;	% 0 for AMDF, 1 for ACF
%		ptOpt=ptOptSet(wObj.fs, wObj.nbits, pfType);
%		ptOpt.mainFun='maxPickingOverPf';
%		ptOpt.usePitchSmooth=0;
%		ptOpt.useVolThreshold=0;
%		ptOpt.useClarityThreshold=0;
%		[pitch, clarity]=pitchTracking(wObj, ptOpt);
%		frameRate=ptOpt.fs/(ptOpt.frameSize-ptOpt.overlap);
%		pitch=[0, 0, pitch];	 % Add leading zeros to start with frame at time zero
%		pitch(pitch==0)=nan;
%		waveCurve4cooledit(wObj, pitch, frameRate);

%		Roger Jang, 20120706

if nargin<1, selfdemo; return; end
if nargin<4, outWaveFile=[tempname, '.wav']; end

% ====== Make "curve" between 0 and 0.8
curveMax=nanmax(curve);
curveMin=nanmin(curve);
curve=0.8*(curve-curveMin)/(curveMax-curveMin);
curveLength=length(curve);
% ====== Creative the interpolated curve
frameTime1=0;								% Time of the first frame
frameTime2=(length(curve)-1)/frameRate;		% Time of the last frame
frameTime=linspace(frameTime1, frameTime2, curveLength)';
sampleTime=(0:frameTime2*wObj.fs)/wObj.fs;
interpolatedCurve=interp1(frameTime, curve, sampleTime);
leng=min(length(wObj.signal), length(interpolatedCurve));
wObj.signal=wObj.signal(1:leng);
interpolatedCurve=interpolatedCurve(1:leng);
wObj2=wObj;
wObj2.signal=[wObj.signal(:), interpolatedCurve(:)];	% Channel 1 is the original signals, channel 2 is the scaled pitch
wObj2file(wObj2, outWaveFile);
dos(['start ', outWaveFile]);		% Start the wav file using cooledit

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
