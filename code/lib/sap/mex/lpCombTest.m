clear all; close all;

fprintf('Compiling lpCombMex.cpp...\n');
mex -I/users/jang/c/lib/audio lpCombMex.cpp /users/jang/c/lib/audio/audioEffect.cpp

fprintf('Testing lpCombMex.mex...\n');
[x, fs, nbits]=wavread('khair.wav');
%x=x(1:3*fs);
prm=lpCombPrmSet;
prm.delay=0.2*fs;

tic
y1=lpComb(x, prm);
time1=toc;
fprintf('time1=%.2f\n', time1);

tic
y2=lpCombMex(x, prm);
time2=toc;
fprintf('time2=%.2f\n', time2);

fprintf('time ratio = %f/%f = %f\n', time1, time2, time1/time2);
fprintf('max(abs(y1-y2) = %f\n', max(abs(y1(:)-y2(:))));
