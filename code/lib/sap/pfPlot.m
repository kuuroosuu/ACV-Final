function pfPlot(pfMat, frameTime, A)
% pfPlot: PF plot, with DP optional path
%	Usage1: pfPlot(pfMat, frameTime, A)
%		where A is obtained from wave2pitchByDp.m
%	Usage2: pfPlot(pfMat, frameTime, pp)
%		where pp is the pitch point on PF matrix.

%	Roger Jang, 20050113, 20100604

if nargin<2, frameTime=1:size(pfMat,2); end
if nargin<3, A=[]; end

imagesc(frameTime, 1:size(pfMat,1), pfMat); shading flat; axis xy;
title('PF matrix');
if isstruct(A)
	for i=1:length(A)
		for j=1:length(A(i).ppc)
			if A(i).ppc(j).index>0
				line(frameTime(i), A(i).ppc(j).index, 'marker', '.', 'color', 'w');
			end
		end
		if isfield(A, 'pp')
			if A(i).pp>0
				line(frameTime(i), A(i).pp, 'marker', 'o', 'color', 'r');
			end
		end
	end
end

if isnumeric(A)		% A is actually pp
	for i=1:length(A)
		line(frameTime(i), A(i), 'marker', '.', 'color', 'w');
	end
end