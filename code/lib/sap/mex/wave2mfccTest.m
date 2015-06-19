fprintf('Compiling wave2mfccMex.cpp...\n');
mex -I/users/jang/c/lib -I/users/jang/c/lib/wave -I/users/jang/c/lib/audio -I/users/jang/c/lib/mfcc -I/users/jang/c/lib/utility wave2mfccMex.cpp \users\jang\c\lib\audio\audio.cpp \users\jang\c\lib\wave\waveRead4pda.cpp \users\jang\c\lib\utility\utility.cpp \users\jang\c\lib\mfcc\CMEL.cpp \users\jang\c\lib\mfcc\CFEABUF.cpp \users\jang\c\lib\mfcc\CSigP.cpp -output wave2mfccMex

addpath /users/jang/matlab/toolbox/sap

waveFile='yankeeDoodle_8k8b.wav';
%waveFile='ivegotrhythm_16k16b.wav';
[y, fs, nbits]=wavReadInt(waveFile);
mfcc=wave2mfccMex(y, fs, nbits, '../mfcc12.cfg');
mesh(mfcc);