% Check one-to-one correspondence between wave files and their "answer.txt".

directory='D:\users\jang\matlab\toolbox\asr\application\anYuan\waveData\簡國慶2';
if isempty(which('findcellstr')); addpath('d:/users/jang/matlab/toolbox/utility'); end
answerFile=[directory, '/answer.txt'];
[songInList, text]=textread(answerFile, '%s %s');
for i=1:length(songInList)
	songInList{i}=[songInList{i}, '.wav'];
end
temp=dir([directory, '/*.wav']);
songInDisk={temp.name};

index=[];
fprintf('In %s, but not in disk:\n', answerFile);
for i=1:length(songInList)
	temp=findcellstr(songInDisk, songInList{i});
	if isempty(temp)
		fprintf('%s in list but not in disk!\n', songInList{i});
		index=[index, i];
	end
end
songInList(index)=[];
text(index)=[];

% 產生新的 answer.txt
%fid=fopen('answer.txt', 'w');
%for i=1:length(songInList)
%	fprintf(fid, '%s\t%s\r\n', songInList{i}(1:end-4), text{i});
%end
%fclose(fid);

index=[];
fprintf('In disk, but not in %s:\n', answerFile);
for i=1:length(songInDisk)
	temp=findcellstr(songInList, songInDisk{i});
	if isempty(temp)
		fprintf('%s in disk but not in list!\n', songInDisk{i});
		index=[index, i];
	end
end