mex waveResampleMex.cpp

addpath d:/users/jang/matlab/toolbox/audio

fileName='主人下馬客在船.wav';
[y, fs1, nbits]=wavReadInt(fileName);
fs2=9131;

% The following commands are synchronous if waveResample.m use the method of "linear" for interpolation.
z1 = waveResample(y, fs1, fs2);
z2 = waveResampleMex(y, fs1, fs2);

time=(1:length(y))/fs1;
time2=(1:length(z2))/fs2;

subplot(2,1,1)
plot(time, y, '.-', time2, z1, '.-', time2, z2, '.-');
grid on
legend('Original', 'reSample.m', 'reSampleMex.dll');
subplot(2,1,2);
plot(1:length(z1), z1-z2, '.-', 1:length(z1), round(z1)-z2, '.-');
grid on