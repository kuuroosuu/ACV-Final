function [insertCount, deleteCount]=insertDeleteCount(onset1, onset2, tolerance, plotOpt)
%insertDeleteCount: Compute insertion and deletion approximately.
%	Usage: [insertCount, deleteCount]=insertDeleteCount(onset1, onset2, tolerance, plotOpt)
%		onset1: system-identified onsets
%		onset2: ground-truth onsets
%
%	For example:
%		onset1=[1.2, 4, 10.2, 15.5, 23];
%		onset2=[1, 4.5, 7, 10, 16, 20];
%		tolerance=1;
%		plotOpt=1;
%		[insertCount, deleteCount]=insertDeleteCount(onset1, onset2, tolerance, plotOpt);
%		fprintf('insertCount=%d, deleteCount=%d\n', insertCount, deleteCount);

%	Roger Jang, 20090608

if nargin<1, selfdemo; return; end
if nargin<4, plotOpt=0; end

onset1=onset1(:);
onset2=onset2(:);
m=length(onset1);
n=length(onset2);

distMat=abs(onset1*ones(1, n)-(onset2*ones(1, m))');
insertCount=sum(min(distMat')>tolerance);
deleteCount=sum(min(distMat)>tolerance);

if plotOpt
	minPos=min(min(onset1), min(onset2));
	maxPos=max(max(onset1), max(onset2));
	plot([minPos, maxPos], [0 0]);
	line(onset1, 0*onset1-0.1, 'marker', '^', 'color', 'k', 'linestyle', 'none');
	line(onset2, 0*onset2+0.1, 'marker', 'v', 'color', 'r', 'linestyle', 'none');
	set(gca, 'ylim', 3*[-1, 1]);
end

% ====== Self demo
function selfdemo
onset1=[1.2, 4, 10.2, 15.5, 23];
onset2=[1, 4.5, 7, 10, 16, 20];
tolerance=1;
plotOpt=1;
[insertCount, deleteCount]=insertDeleteCount(onset1, onset2, tolerance, plotOpt);
fprintf('insertCount=%d, deleteCount=%d\n', insertCount, deleteCount);