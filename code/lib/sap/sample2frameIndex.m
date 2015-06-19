function frameIndex=sample2frameIndex(sampleIndex, frameSize, overlap)
%sample2frameIndex: Sample to frame index conversion
%	Usage: frameIndex=sample2frameIndex(sampleIndex, frameSize, overlap)

%	Roger Jang, 20041016

if nargin<1; selfdemo; return; end

frameIndex=1+round((sampleIndex-1)/(frameSize-overlap));

% ====== Self demo
function selfdemo
sampleIndex=[3 5 7];
frameSize=5;
overlap=3;
fprintf('sampleIndex = %s\n', mat2str(sampleIndex));
fprintf('frameSize = %d\n', frameSize);
fprintf('overlap = %d\n', overlap);
fprintf('After calling "frameIndex=sample2frameIndex(sampleIndex, frameSize, overlap)":\n');
frameIndex=sample2frameIndex(sampleIndex, frameSize, overlap);
fprintf('frameIndex = %s\n', mat2str(frameIndex));
