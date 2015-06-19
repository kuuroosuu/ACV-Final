function volTh=getVolTh(volume, frameSize, method)
% getVolTh: Get volume threshold
%	Usage: volTh=getVolTh(volume, method)

%	Roger Jang, 20050105

if nargin<2, frameSize=256; end
if nargin<3, method=1; end

switch method
	case 1
		sortedVol=sort(volume);
		volTh=sortedVol(round(length(volume)/4));
	case 2
		testNum=5;
		volTh = 4*mean(volume(1:testNum));
		if (volTh<frameSize/2) | (volTh>10*frameSize)
			volTh = 2*frameSize;
		end
	case 3
		volTh=10*min(volume);
	case 4
		volTh=max(volume)/10;
	otherwise
		error('Unknown method!');
end