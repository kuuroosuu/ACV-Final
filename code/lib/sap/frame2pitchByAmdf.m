function [pitch, clarity, pitchIndex]=frame2pitchByAmdf(frame, ptOpt, plotOpt)
% frame2pitchByAmdf: Frame to pitch conversion using AMDF
%	Usage: [pitch, clarity, pitchIndex]=frame2pitchByAmdf(frame, ptOpt, plotOpt)
%
%	See also FRAME2ACF, FRAME2AMDF.

%	Roger Jang 20070209

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

amdf=frame2amdfMex(frame, length(frame), 2);	% This is the major difference from frame2pitchByAcf!
[pitch, clarity, pitchIndex]=pf2pitch(amdf, ptOpt);

if plotOpt
	subplot(2,1,1); plot(frame); axis tight; title('Frame');
	subplot(2,1,2); [pitch, clarity, pitchIndex]=pf2pitch(amdf, ptOpt, plotOpt); title('AMDF');
end

% ====== Self demo
function selfdemo
waveFile='greenOil.wav';
[y, fs, nbits]=wavReadInt(waveFile);
framedY=buffer2(y, 256, 0);
frame=framedY(:, 250);
ptOpt=ptOptSet(fs, nbits, 0);
plotOpt=1;
[pitch, clarity, pitchIndex]=frame2pitchByAmdf(frame, ptOpt, plotOpt);