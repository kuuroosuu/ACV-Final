function distance=ioiDistance(queryIoi, dbIoi)
% qbtDistance: Compute the distance between query and database IOIs
%	Usage: ioiDistance(queryIoi, dbIoi)

if length(queryIoi)<2, distance=inf; return; end
if length(dbIoi)<2, distance=inf; return; end
queryIoi=queryIoi(:);
dbIoi=dbIoi(:);
queryIoi(queryIoi==0)=[];
dbIoi(dbIoi==0)=[];
queryIoiRatio=queryIoi(2:end)./queryIoi(1:end-1);
dbIoiRatio=dbIoi(2:end)./dbIoi(1:end-1);
noteNum=min(length(queryIoiRatio), length(dbIoiRatio));
distance=sum(abs(queryIoiRatio(1:noteNum)-dbIoiRatio(1:noteNum)))/noteNum;