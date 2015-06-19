function lineStyle=getLineStyle(index)
% getLineStyle: Get line style from a rotating palette
%	Usage: lineStyle=getLineStyle(index)

% Roger Jang, 20071020

if nargin<1, index=1; end

allLineStyle={'-', '--', ':', '-.'};
lineStyle=allLineStyle{mod(index-1, length(allLineStyle))+1};