function [cueTime, msg] = wavreadCue(wavFile, plotOpt)
% [cueTime, msg] = wavreadCue(wavFile, plotOpt)
% cueTime: time sequence
% msg : message 
% wavFile : wavefile
% Hong-Ru Lee, 20081027

if nargin<1, selfdemo; return; end
if nargin<2, plotOpt=0; end

% ====== Get fs
[y, fs]=wavread(wavFile);

% ====== Get file contents
fid=fopen(wavFile);
contents=fread(fid)';
fclose(fid);

% ====== find Cue Chunk ID
index = findstr(contents, 'cue ');

% ====== Number of Cue Mark
tapNum = contents(index+8);

for i=1:tapNum
	tmp=[];
    i;
	for j=4:-1:1
		if length(num2str(dec2hex(contents(index+16+24*(i-1)+j-1))))==1
			tmp = [tmp '0' dec2hex(contents(index+16+24*(i-1)+j-1))];
		else
			tmp = [tmp dec2hex(contents(index+16+24*(i-1)+j-1))];
		end
   	end
	tapDuration(i) = hex2dec(tmp);
end

if isempty(index),
	cueTime=[];
	msg='Not a cue chunk of wave';
else
	cueTime=sort(tapDuration)/fs;
	msg='correct';
end

if plotOpt
	time=(1:length(y))/fs;
	plot(time, y);
	for i=1:length(cueTime)
		pos=cueTime(i);
		line(pos*[1 1], [-1 1], 'color', 'r');
	end
end

% ====== Self demo
function selfdemo
wavFile='0107.wav';
plotOpt=1;
wavReadCue(wavFile, plotOpt);