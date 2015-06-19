function y2=waveResample(y, fs, fs2, plotOpt)
%RESAMPLE Resampling of signals using spline interpolation
%	Usage: y2=waveResample(y, fs, fs2, plotOpt)

%	Roger Jang, 20030624

if nargin<1, selfdemo; return; end
if nargin<2, fs=16000; end
if nargin<3, fs2=8000; end	% 預設是降頻到 8KHz
if nargin<4, plotOpt=0; end

if fs==fs2
	y2=y;
	return;
end

time=(0:length(y)-1)/fs;
time2=(0:floor((length(y)-1)*fs2/fs))/fs2;
y2=interp1(time, y, time2, 'spline');

if plotOpt
	plot(time, y, '.-', time2, y2, '.-');
	xlabel('Time (second)');
	title(sprintf('Original signal y (fs=%d) and resampled signal y2 (fs=%d)', fs, fs2));
	legend(sprintf('Original signal y (fs=%d)', fs), sprintf('Resampled signal y2 (fs=%d)', fs2));
end

% ====== selfdemo
function selfdemo
fileName='主人下馬客在船.wav';
[y, fs, nbits]=wavReadInt(fileName);
fs2=8000;
y2=feval(mfilename, y, fs, fs2, 1);

fprintf(sprintf('Press any key to hear the original sound (fs=%d):', fs)); pause; fprintf('\n');
sound(y/(2^nbits/2), fs);
fprintf(sprintf('Press any key to hear the resampled sound (fs=%d):', fs2)); pause; fprintf('\n');
sound(y2/(2^nbits/2), fs2);