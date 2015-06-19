close all;
clear all
fprintf('Compiling frame2acfMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility -I/users/jang/c/lib/wave frame2acfMex.cpp \users\jang\c\lib\audio\audio.cpp \users\jang\c\lib\utility\utility.cpp \users\jang\c\lib\wave\waveRead4pda.cpp
fprintf('Compiling frame2acfIntMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility -I/users/jang/c/lib/wave frame2acfIntMex.cpp \users\jang\c\lib\audio\audio.cpp \users\jang\c\lib\utility\utility.cpp \users\jang\c\lib\wave\waveRead4pda.cpp
fprintf('Finish compiling\n');
dos('copy /y frame2acf*.mex* ..\private');
addpath ..

[y, fs, nbits]=wavReadInt('七月七日長生殿.wav');
%[y, fs, nbits]=wavread('七月七日長生殿.wav');
ptOpt=ptOptSet(fs, nbits);
frameSize=ptOpt.frameSize;
overlap=ptOpt.overlap;
frameMat=enframe(y, frameSize, overlap);
frame=frameMat(:, 90);
subplot(3,1,1);
plot(frame);
subplot(3,1,2);
acf1=frame2acf(frame);
acf2=frame2acfMex(frame);
acf3=frame2acfIntMex(frame);
plot([acf1 acf2 acf3]);
legend('acf1', 'acf2', 'acf3');
subplot(3,1,3);
plot([acf1-acf2, acf2-acf3]);
 legend('acf1-acf2', 'acf2-acf3');
diff12=sum(abs(acf1-acf2));
diff23=sum(abs(acf2-acf3));
fprintf('Diff between frame2acf() and frame2acfMex() = %f\n', diff12);
fprintf('Diff between frame2acfMex() and frame2acfIntMex() = %f\n', diff23);
