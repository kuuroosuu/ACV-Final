function threshold=getVolumeThreshold(volume, testNum, divisor)
% Get the threshold of volume for rejecting a pitch

if nargin<2, testNum=5; end
if nargin<3, divisor=4; end

threshold=floor((sum(volume(1:testNum))+floor(divisor/2))/divisor);

if (threshold<256/2) | (threshold>10*256)
	threshold=2*256;
end