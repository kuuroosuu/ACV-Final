function cHitTime=hitTimeDetect(wObj, opt, showPlot)
% hitTimeDetect: Detect the hit time of a give wave object
%
%	Usage:
%		cHitTime=hitTimeDetect(wObj, opt, showPlot)
%
%	Example:
%		waveFile='song01s5.wav';
%		wObj=waveFile2obj(waveFile);
%		opt=hitTimeDetect('defaultOpt');
%		cHitTime=hitTimeDetect(wObj, opt, 1);
%		tempWaveFile=[tempname, '.wav'];
%		cueLabelAdd(waveFile, tempWaveFile, cHitTime);
%		dos(['start ', tempWaveFile]);

%	Roger Jang, 20130707

if nargin<1, selfdemo; return; end
% ====== Set the default options
if nargin==1 && ischar(wObj) && strcmpi(wObj, 'defaultOpt')
	cHitTime.oscOpt=wave2osc('defaultOpt');
	cHitTime.hitTimeOpt=osc2hitTime('defaultOpt');
	return
end
if nargin<2|isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, showPlot=0; end

if ischar(wObj), wObj=waveFile2obj(wObj); end
wObj.signal=mean(wObj.signal, 2);			% Converted into mono

osc=wave2osc(wObj, opt.oscOpt, showPlot);
if showPlot, figure; end
cHitTime=osc2hitTime(osc, opt.hitTimeOpt, showPlot);
%cueLabelAdd(waveFile, waveData(i).waveFileWithCues, cHitTime);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
