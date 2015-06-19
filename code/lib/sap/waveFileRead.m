function [y, fs, nbits]=waveFileRead(waveFile, fs2, nbits2)
% waveFileRead: Read a wave file and convert the signals to the given format
%	Usage: [y, fs, nbits]=waveFileRead(waveFile, fs2, nbits2)
%	waveFile: wave file
%	fs2: New sample rate
%	nbits2: New bit resolution, either 8 or 16
%	y: Unsigned integers between 0 and 2^nbits2-1
%	fs: Original sample rate
%	nbits: Original bit resolution

%	Roger Jang, 20040103

if nargin==0; selfdemo; return; end

[y, fs, nbits]=wavread(waveFile);
y=y*(2^nbits/2)+2^nbits/2;		% 轉換成 unsigned integer (between 0 and 2^nbits-1)
if (nargin>=2) & (fs2~=fs)
	fprintf('Resample of "%s" from %g to %g...\n', waveFile, fs, fs2);
	y=floor(resample(y, fs, fs2)+0.5);	% Resample to fs2
end
if (nargin>=3) & (nbits2~=nbits)
	fprintf('Change bit resolution of "%s" from 16-bit to 8-bit...\n', waveFile);
	if (nbits==16) & (nbits2==8)
		y=floor((y+2^7)/2^8);	% Change bit resolution from 16 to 8 bits
	end
	if (nbits==8) & (nbits2==16)
		y=y*2^8;		% Change bit resolution from 8 to 16 bits
	end
end


function selfdemo
waveFile='主人下馬客在船.wav';
[y1, fs, nbits]=feval(mfilename, waveFile);
t1=(0:length(y1)-1)/fs;
fs2=8000;
[y2, fs, nbits]=feval(mfilename, waveFile, fs2);
t2=(0:length(y2)-1)/fs2;
plot(t1, y1, '.-', t2, y2, 'o-');
legend(['fs=', int2str(fs)], ['fs=', int2str(fs2)]);