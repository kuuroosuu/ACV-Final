function gtBeat=btGtRead(gtBeatFile)
% btGtRead: Beat tracking ground-truth reading

[parentDir, mainName, extName]=fileparts(gtBeatFile);
if ~strcmp(extName, 'txt')
	gtBeatFile=strrep(gtBeatFile, extName, '.txt');
end

if ~exist(gtBeatFile, 'file')
	gtBeat=[];
	return
end

gtBeat{1}=asciiRead(gtBeatFile);
return

contents=textread(gtBeatFile, '%s', 'delimiter', '\n', 'bufsize', 20000);
for j=1:length(contents)
	temp=split(contents{j}, 9);
	if isempty(temp{end}), temp(end)=[]; end	% Delete the last element if it's empty
	gtBeat{j}=str2double(temp);
end