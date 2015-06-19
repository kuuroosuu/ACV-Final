function y2=resample(y, fs, fs2)
%RESAMPLE Resampling of signals
%	Usage: y2=resample(y, fs, fs2)

%	Roger Jang, 20020721

if nargin==0, selfdemo; return; end
if nargin<3, fs2=8000; end	% 預設是降頻到 8KHz

time=(1:length(y))/fs;
time2=(1:floor(length(y)*fs2/fs))/fs2;
y2=interp1(time, y, time2);

% ====== selfdemo
function selfdemo
fileName='water.wav';
[y, fs]=wavread(fileName);
fs2=8000;
y2=feval(mfilename, y, fs, fs2);

subplot(2,1,1);
time=(1:length(y))/fs;
time2=(1:length(y2))/fs2;

subplot(2,1,1);
plot(time, y, time, y, 'r.');
title(sprintf('Original signal y (fs=%d)', fs));
axis([-inf, inf, -1, 1]);
subplot(2,1,2);
plot(time2, y2, time2, y2, 'r.');
title(sprintf('Resampled signal y2 (fs2=%d)', fs2));
axis([-inf, inf, -1, 1]);
xlabel('Time (second)');

fprintf('Press any key to hear the original sound:'); pause; fprintf('\n');
sound(y, fs);
fprintf('Press any key to hear the resampled sound:'); pause; fprintf('\n');
sound(y2, fs2);