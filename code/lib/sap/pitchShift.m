function wObj2=pitchShift(wObj, opt, plotOpt)
% pitchShift: Pitch Shifting
%
%	Usage:
%		wObj2=pitchShift(wObj)
%		wObj2=pitchShift(wObj, opt)
%		wObj2=pitchShift(wObj, opt, plotOpt)
%
%	Description:
%		wObj2=pitchShift(wObj, opt) return the wave object after pitch shifting.
%			wObj: Input wave object
%			opt: Options for pitch shifting
%				opt.pitchShiftAmount: Amount of pitch shift (in semitone)
%				opt.method: 'wsola' or 'phaseVocoder'
%			wObj2: Output wave object
%
%	Example:
%		% === Pitch shift using wsola
%		wObj=waveFile2obj('what_movies_have_you_seen_recently.wav');
%		opt=pitchShift('defaultOpt');
%		opt.method='wsola';
%		wObj2=pitchShift(wObj, opt, 1);
%		% === Pitch shift using phase vocoder
%		wObj=waveFile2obj('what_movies_have_you_seen_recently.wav');
%		opt=pitchShift('defaultOpt');
%		opt.method='phaseVocoder';
%		figure;
%		wObj2=pitchShift(wObj, opt, 1);
%
%	See also wsola.

%	Roger Jang, 20120607, 20120911

if nargin<1, selfdemo; return; end
if nargin==1 && ischar(wObj) && strcmpi(wObj, 'defaultOpt')	% Set the default options
	wObj2.pitchShiftAmount=-6;	% Amount for pitch shift (in semitones)
	wObj2.method='wsola';		% Method for pitch shifting
	return
end
if nargin<2||isempty(opt), opt=feval(mfilename, 'defaultOpt'); end
if nargin<3, plotOpt=0; end
if ischar(wObj), wObj=waveFile2obj(wObj); end	% wObj is actually the wave file name

if opt.pitchShiftAmount==0, wObj2=wObj; return; end
if size(wObj.signal, 2)>1, error('Stereo sound is not supported in %s!\n', mfilename); end

durationRatio=2^(opt.pitchShiftAmount/12);
switch(opt.method)
	case 'wsola'
		wsolaOpt=wsola('defaultOpt');
		wsolaOpt.durationRatio=durationRatio;
		wObj2=wsola(wObj, wsolaOpt);
	case 'phaseVocoder'
		wObj2=wObj;
		wObj2.signal=pvoc(wObj2.signal, 1/durationRatio, 1024);
	otherwise
		error('Unknown method=%s in %s!\n', opt.method, mfilename);
end
[n, d]=rat(durationRatio);
wObj2.signal=resample(wObj2.signal, d, n);

if plotOpt
	maxTime=max(length(wObj.signal)/wObj.fs, length(wObj2.signal)/wObj2.fs);
	subplot(2,1,1); plot((1:length(wObj.signal))/wObj.fs, wObj.signal);
	title('Original');
	axis([0, maxTime, -1, 1]);
	wavePlayButton(wObj);
	subplot(2,1,2); plot((1:length(wObj2.signal))/wObj2.fs, wObj2.signal);
	title(sprintf('Synthesized (method=%s, pitchShiftAmount=%g)', opt.method, opt.pitchShiftAmount));
	axis([0, maxTime, -1, 1]);
	xlabel('Time (sec)');
	wavePlayButton(wObj2);
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
