function endPointList=endPointRead(endPointFile)
% endPointRead: Read a given endpoint file
%	Usage: endPointRead(endPointFile)

%	Roger Jang, 20031120

if nargin<1, endPointFile='waveFile\endpoint.txt'; end

contents = textread(endPointFile,'%s','delimiter','\n','whitespace','');

fileNum=length(contents);
for i=1:fileNum
	items=split(contents{i}, [9, 32]);
	if length(items)==0	% ¸õ±¼ªÅ¦æ
		continue;
	end
	endPointList(i).fileName=items{1};
	endPointList(i).startIndex=eval(items{2});
	endPointList(i).endIndex=eval(items{3});
end