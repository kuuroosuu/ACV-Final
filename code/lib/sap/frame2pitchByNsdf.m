function [pitch, clarity, pitchIndex]=frame2pitchByNsdf(frame, ptOpt, plotOpt)
% frame2pitchByNsdf: Frame to pitch conversion using NSDF
%	Usage: [pitch, clarity, pitchIndex]=frame2pitchByNsdf(frame, ptOpt, plotOpt)
%
%	See also FRAME2ACF, FRAME2SMDF.

%	Roger Jang 20070209

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

nsdf=frame2nsdfMex(frame);
[pitch, clarity, pitchIndex]=nsdf2pitch(nsdf, ptOpt);

if plotOpt
	subplot(2,1,1); plot(frame); axis tight; title('Frame');
	subplot(2,1,2); [pitch, clarity, pitchIndex]=nsdf2pitch(nsdf, ptOpt, plotOpt); title('NSDF');
end

% ====== Self demo
function selfdemo
waveFile='greenOil.wav';
[y, fs, nbits]=wavReadInt(waveFile);
framedY=buffer2(y, 256, 0);
frame=framedY(:, 250);
ptOpt=ptOptSet(fs, nbits);
plotOpt=1;
[pitch, clarity, pitchIndex]=frame2pitchByNsdf(frame, ptOpt, plotOpt);