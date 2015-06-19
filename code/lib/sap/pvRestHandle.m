function outputPv=pvRestHandle(inputPv, useRest, plotOpt)
% pvRestHandle: Handle rest in PV for QBSH
%	Usage: outputPv=pvRestHandle(inputPv, useRest)
%		inputPv: input PV (pitch vector)
%		useRest: 0 ===> remove rests, 1 ===> replace a rest by its previous non-rest pitch

%	Roger Jang, 20040613

if nargin<1, selfdemo; return; end
if nargin<3, plotOpt=0; end

outputPv=inputPv;

% ====== Cut out leading/trailing rests
while outputPv(1)==0; outputPv(1)=[]; end
while outputPv(end)==0; outputPv(end)=[]; end

% ====== Handle in-between rests
if ~useRest
	outputPv(outputPv==0)=[];		% 砍掉所有零點
else
	for j=2:length(outputPv)	% 碰到零點，延長前一個音
		if outputPv(j-1)~=0 & outputPv(j)==0
			outputPv(j)=outputPv(j-1);
		end
	end
end

if plotOpt==1
	plot(outputPv, '.-');
end

% ====== Self demo
function selfdemo
load london_bridge.pv
inputPv=london_bridge;
tempPv=inputPv
tempPv(tempPv==0)=nan;
subplot(3,1,1);
plot(tempPv, '.-'); title('Original PV');
subplot(3,1,2);
outputPv=pvRestHandle(inputPv, 0, 1); title('useRest=1');
subplot(3,1,3);
outputPv=pvRestHandle(inputPv, 1, 1); title('useRest=0');
xlabel('Frame index');