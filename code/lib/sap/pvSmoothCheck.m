function pvSmoothCheck(pvDir, opt)
% pvSmoothCheck: Identify possible pitch errors within a folder of pv files.
%
%	Usage:
%		pvSmoothCheck(pvDir, opt)
%
%	Description:
%		The errors identified by the program include:
%			errorCount(1): No. of occurrences of “0*0” patterns, where * indicates a nonzero pitch.
%			errorCount(2): No. of occurrences of “*0*” patterns, where * indicates a nonzero pitch.
%			errorCount(3): No. of occurrences of two neighboring pitch points with difference higher than 5 semitones.
%			errorCount(4): No. of occurrences of pitch points higher than 70 (for speech) or 80 (for singing) semitones.
%			errorCount(5): No. of pitch segments minus the number of syllables. (If this is negative, set it to zero.) (Only for speech)
%		Use this program to check your PV files. Make sure your PV files do not suffer from any of the above errors.
%
%	Example:
%		pvDir='D:\dataSet\tangPoem\2011-音訊處理與辨識';			% Speech
%		pvDir='D:\dataSet\childSong\waveFile\2011-音訊處理與辨識';		% Singing
%		opt.type='singing';
%		opt.maxPitch=80;
%		opt.maxPitchDiff=5;
%		opt.outputDir=tempname;
%		pvSmoothCheck(pvDir, opt);

%	Roger Jang, 20101114, 20110503

if nargin<1, selfdemo; return; end
if nargin<2 || isempty(opt)
	opt.type='singing';
	opt.maxPitch=80;
	opt.maxPitchDiff=5;
	opt.outputDir=tempname;
end

if exist(opt.outputDir, 'dir')~=7
	fprintf('Creating %s directory...\n', opt.outputDir);
	mkdir(opt.outputDir);
end

pvData=recursiveFileList(pvDir, 'pv');
if isempty(pvData)
	fprintf('No PV files found within %s!\n', pvDir);
	fprintf('Exiting...\n');
	return
end

fprintf('Reading pv and do pv segmentation within %d pv files...\n', length(pvData));
tic
for i=1:length(pvData)
	pvData(i).pv=asciiRead(pvData(i).path);
	pvData(i).segment=segmentFind(pvData(i).pv);
	pvData(i).errorCount=[];
end
fprintf('\tElapsed time = %.2f sec\n', toc);

errorDoc={};

fprintf('Counting the occurrence of 0*0 within %d pv files...\n', length(pvData));
errorDoc{end+1}='No. of occurrences of 0*0';
tic
for i=1:length(pvData)
	pv=pvData(i).pv;
	pvData(i).errorCount(end+1)=sum(~pv(1:end-2) & pv(2:end-1) & ~pv(3:end));
end
fprintf('\tElapsed time = %.2f sec\n', toc);

if strcmp(opt.type, 'speech')
fprintf('Counting the occurrence of *0* within %d pv files...\n', length(pvData));
errorDoc{end+1}='No. of occurrences of *0*';
tic
for i=1:length(pvData)
	pv=pvData(i).pv;
	pvData(i).errorCount(end+1)=sum(pv(1:end-2) & ~pv(2:end-1) & pv(3:end));
end
fprintf('\tElapsed time = %.2f sec\n', toc);
end

fprintf('Counting the occurrence of highly different neighboring pitch within %d pv files...\n', length(pvData));
errorDoc{end+1}=sprintf('No. of occurrences of neighboring pitch difference larger than %d', opt.maxPitchDiff);
tic
for i=1:length(pvData)
	pv=pvData(i).pv;
	segment=pvData(i).segment;
	errorCount=0;
	for j=1:length(segment)
		theCount=sum(abs(diff(pv(segment(j).begin:segment(j).end)))>=opt.maxPitchDiff);
		errorCount=errorCount+theCount;
	end
	pvData(i).errorCount(end+1)=errorCount;
end
fprintf('\tElapsed time = %.2f sec\n', toc);

fprintf('Counting the occurrence of excessive high pitch within %d pv files...\n', length(pvData));
errorDoc{end+1}=sprintf('No. of pitch points over %d semitone', opt.maxPitch);
tic
for i=1:length(pvData)
	pvData(i).errorCount(end+1)=sum(pvData(i).pv>=opt.maxPitch);
end
fprintf('\tElapsed time = %.2f sec\n', toc);

if strcmp(opt.type, 'speech')
fprintf('Counting the difference between number of pitch segments and number of syllables within %d pv files...\n', length(pvData));
errorDoc{end+1}='Difference between number of pitch segments and number of syllables.';
tic
for i=1:length(pvData)
	pvData(i).syllableCount=length(strrep(pvData(i).name, '.pv', ''));
	pvData(i).segmentCount=length(pvData(i).segment);
	pvData(i).errorCount(end+1)=max(pvData(i).segmentCount-pvData(i).syllableCount, 0);
end
fprintf('\tElapsed time = %.2f sec\n', toc);
end

for i=1:length(pvData)
	pvData(i).allErrorCount=sum(pvData(i).errorCount);
end

pvData([pvData.allErrorCount]==0)=[];		% Get rid of entry with zero error count
[sorted, index]=sort([pvData.allErrorCount], 'descend');
pvData=pvData(index);
for i=1:length(pvData)
	if pvData(i).errorCount==0, break; end
	pvFile=pvData(i).path;
	fprintf('%d/%d: Error counts = %s, file = %s\n', i, length(pvData), mat2str(pvData(i).errorCount), pvFile);
	pv=asciiRead(pvFile);
	plot(pv, '.-'); title(strPurify(sprintf('PV file=%s\n', pvFile)));
	pvData(i).image=sprintf('%s_%s', pvData(i).parentDir, strrep(pvData(i).name, 'pv', 'png'));
	imagePath=sprintf('%s/%s', opt.outputDir, pvData(i).image);
	if ~exist(imagePath, 'file')
		cmd=sprintf('print(''-dpng'',  ''%s'');', strrep(imagePath, '''', ''''''));
		fprintf('\tCreating %s...\n', imagePath);
		eval(cmd);
	end
%	fprintf('Hit any key to continue...\n'); pause; fprintf('\n');
end
titleStr='Table of error count for pv files';
titleStr=sprintf('%s<ol>', titleStr);
for i=1:length(errorDoc)
	titleStr=sprintf('%s<li>%s', titleStr, errorDoc{i});
end
titleStr=sprintf('%s</ol>', titleStr);
structDispInHtml(pvData, titleStr, {'parentDir', 'path', 'errorCount', 'image'}, {'Owner', 'PV file', 'Error count', 'PV curve'}, {'image'}, [opt.outputDir, '/index.htm']);

% ====== Self demo
function selfdemo
mObj=mFileParse(which(mfilename));
strEval(mObj.example);
