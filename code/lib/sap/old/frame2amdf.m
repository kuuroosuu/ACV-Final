function amdf=frame2amdf(frame, maxShift, plotOpt, mode)
% FRAME2AMDF Frame to AMDF (average of magnitude difference function) for pitch tracking
%	Usage: amdf=amdf(frame, maxShift, plotOpt);
%	amdf: Returned amdf vector

%	Roger Jang 20020404

if nargin<1, selfdemo; return; end
if nargin<2, maxShift=floor(length(frame)/2); end
if nargin<3, plotOpt=0; end

frameSize=length(frame);
amdf=zeros(maxShift, 1);
for i = 2:maxShift
	amdf(i) = sum(abs(frame(1:frameSize-maxShift)-frame(i:frameSize-maxShift+i-1)));	% 平移 frame(1:frameSize-maxShift)
%	amdf(i) = sum(abs(frame(1:maxShift) - frame(i:maxShift+i-1)));				% 平移 frame(1:maxShift)
%	amdf(i) = floor(amdf(i)/(frameSize- maxShift));	% 求平均
end

amdf = round(amdf/(frameSize - maxShift));

if plotOpt
	subplot(2,1,1);
	plot(frame, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('Frame');
	subplot(2,1,2);
	plot(amdf, '.-');
	set(gca, 'xlim', [-inf inf]);
	title('AMDF vector');
end


% ====== Self demo
function selfdemo
waveFile='greenOil.wav';
[y, fs, nbits]=waveFileRead(waveFile);
framedY=buffer2(y, 256, 0);
frame=framedY(:, 250);
plotOpt=1;
feval(mfilename, frame, 128, plotOpt);