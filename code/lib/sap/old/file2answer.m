function answer=file2answer(file, answerFile)
% file2answer: 給一個檔案（可使用絕對路徑），根據同目錄下的 answer.txt 回傳相關的答案
%	Usage:	answer=file2answer(file, answerFile)
%		answer=file2answer(file)

if nargin==0; selfdemo; return; end
if isempty(which('findcellstr')); addpath('d:/users/jang/matlab/toolbox/utility'); end
[parentDir, mainName, extName, version]=fileparts(file);
if isempty(parentDir), parentDir='.'; end
if nargin<2; answerFile=[parentDir, '\answer.txt']; end

[fileNames, answers]=textread(answerFile, '%s\t%s');
index=findcellstr(fileNames, mainName);
if ~isempty(index)
	if length(index)>1, error('More than one answer is returned!'); end
	answer=answers{index};
else
	answer=[];
end

% ===== self demo
function selfdemo
file='D:\users\jang\matlab\toolbox\asr\application\RMC\ivr辨識音檔2\alexxx1.wav';
output=feval(mfilename, file);
fprintf('Input: %s\n', file);
fprintf('Output: %s\n', output);