function freq=pevRead(pevFile)
% rawRead: Read a .pev file from Keele pitch data set
%	Usage: freq=pevRead(pevFile)
%	
%	For example:
%	
%	freq=pevRead('D:\dataSet\KeelePitch\Speech\f1nw0000.pev');

%	Roger Jang, 20051202

contents = textread(pevFile, '%s', 'delimiter', '\n', 'whitespace', '');

freq=[];
for i=14:length(contents)
	if strcmp(contents{i}, 'ELF:')
		break;
	end
	freq=[freq; eval(contents{i})];
end