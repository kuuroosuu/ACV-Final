clear all; close all;

fprintf('Compiling tapDelayLineMex.cpp...\n');
mex -I/users/jang/c/lib/audio tapDelayLineMex.cpp /users/jang/c/lib/audio/audioEffect.cpp

fprintf('Testing tapDelayLineMex.mex...\n');
%x=zeros(100,1); x(1)=1;		% unit impulse signal of length 100
%prm.delay=[10 24 37];
%prm.gain=[0.8 0.6 0.3];
%y=tapDelayLineMex(x, prm);

[x, fs, nbits]=wavread('khair.wav');
%x=x(1:3*fs);
prm=tapDelayLinePrmSet;
prm.delay=[0.1, 0.3, 0.5]*fs;
prm.gain=[0.8, 0.5, 0.2];

tic
y1=tapDelayLine(x, prm);
time1=toc;
fprintf('time1=%.2f\n', time1);

tic
y2=tapDelayLineMex(x, prm);
time2=toc;
fprintf('time2=%.2f\n', time2);

fprintf('time ratio = %f/%f = %f\n', time1, time2, time1/time2);
fprintf('max(abs(y1-y2) = %f\n', max(abs(y1(:)-y2(:))));
