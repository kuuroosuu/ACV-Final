function [ticks, tick]=tickSequence(opt, showPlot)
%tickSequence: Create a tick sequence for a metronome
%
%	Usage:
%		ticks=tickSequence(opt, showPlot)
%
%	Example:
%		opt=tickSequence('defaultOpt');
%		ticks=tickSequence(opt, 1);
%		sound(ticks, opt.fs);
%		%fprintf('Saving to tickSequence.wav...\n');
%		%wavwrite(ticks, opt.fs, 'tickSequence');

%	Roger Jang, 20130429

if nargin<1, selfdemo; return; end
% ====== Set the default options
if nargin==1 && ischar(opt) && strcmpi(opt, 'defaultOpt')
	ticks.bpm=120;
	ticks.beatPerMeasure=4;
	ticks.freq01=600;
	ticks.freq02=800;
	ticks.duration=10;
	ticks.fs=16000;
	return
end
if nargin<2, showPlot=0; end

tickSingleOpt=tickSingle('defaultOpt');
tickSingleOpt.freq=opt.freq01;
tickSingleOpt.fs=opt.fs;
tick01=tickSingle(tickSingleOpt);
tickSingleOpt.freq=opt.freq02;
tick02=tickSingle(tickSingleOpt);

ibi=round(60/opt.bpm*opt.fs);			% ibi=inter-beat interval (in sample point)
if ibi<0.1*length(tick01), error('IBI too short!'); end
aPeriod=zeros(ibi, 1);
periodNum=ceil(opt.duration*opt.fs/ibi);
mat=[];
for i=1:periodNum
	if mod(i, opt.beatPerMeasure)==1
		aPeriod(1:length(tick02))=tick02';
	else
		aPeriod(1:length(tick01))=tick01';
	end
	mat=[mat, aPeriod];
end
%mat=repmat(aPeriod, 1, periodNum);
ticks=mat(:)';
ticks=ticks(1:opt.duration*opt.fs);

if showPlot
	time=(0:length(tick01)-1)/opt.fs;
	subplot(221); plot(time, tick01); title('Low tick'); xlabel('Time (sec)');
	subplot(222); plot(time, tick02); title('High tick'); xlabel('Time (sec)');
	subplot(212);
	time=(0:opt.duration*opt.fs-1)/opt.fs;
	plot(time, ticks); title('Tick sequence'); xlabel('Time (sec)');
end

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
