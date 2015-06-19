function [y, fs, nbits]=mp3read(mp3file, playOpt)
% mp3read: Read mp3 to int16 output
%	Usage: [y, fs, nbits]=mp3read(mp3file, playOpt)

if nargin<1, selfdemo; return; end
if nargin<2, playOpt=0; end

fs=44100;
nbits=16;
% Convert mp3 file to raw format
cmd=['d:\users\jang\matlab\toolbox\sap\mp32wav "', mp3file, '"'];
fprintf('Executing ¡u%s¡v...\n', cmd);
dos(cmd);
fileName='temp.raw';	% Raw data file

% Make sure the file conversion is complete
%for i=1:1000,
%	fid=fopen(fileName, 'w');
%	if fid>0,
%		break;
%	end
%	fprintf('file not ready\n');
%	pause(1);
%	fclose(fid);
%end

% Read temp.raw
fid=fopen(fileName, 'r');
fprintf('Reading data from "%s"...\n', fileName);
s=fread(fid, inf, '*short');
fclose(fid);

%fprintf('Reshaping...\n');
y=reshape(s, 2, length(s)/2)';

%fprintf('Converting to a double between -1 and 1...\n');
%y=double(y)/(2^15);

fprintf('Deleting "temp.raw"...\n');
delete temp.raw

if playOpt
	duration=10;
	fprintf('Hit return to play 10 seconds of "%s"...', mp3file); pause
	sound(double(y(1:duration*fs))/(2^15), fs);
end

% ====== Self demo
function selfdemo
file='Charlene - I''ve Never Been To Me.mp3';
[y, fs, nbits]=feval(mfilename, file, 1);
