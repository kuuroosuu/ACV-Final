% ====== 直接加入 Utility Toolbox 和 SAP Toolbox
%addpath d:/users/jang/matlab/toolbox/utility
%addpath d:/users/jang/matlab/toolbox/sap

% ====== 經由對話窗加入 Utility 和 SAP Toolbox
if exist('recursiveFileList')==0
	fprintf('請先下載 Utility Toolbox（http://mirlab.org/jang/matlab/toolbox/utility.rar）並加入搜尋路徑！\n');
	dirName=uigetdir('', '請點選 Utility Toolbox 之目錄');
	addpath(dirName);
end
if exist('frameFlip')==0
	fprintf('請先下載 SAP Toolbox（http://mirlab.org/jang/matlab/toolbox/sap.rar）並加入搜尋路徑！\n');
	dirName=uigetdir('', '請點選 SAP Toolbox 之目錄');
	addpath(dirName);
end

% ====== 檢查 MATLAB 版本
%desiredVersion='6.5';
%matlabVersion=version;
%if ~strcmp(matlabVersion(1:3), desiredVersion)
%	fprintf('建議使用版本：%s\n', desiredVersion);
%end