close all; clear all
addMyPath;
fprintf('Compiling pitch2waveMex.cpp...\n');
mex pitch2waveMex.cpp \users\jang\c\lib\audio\audio.cpp \users\jang\c\lib\utility\utility.cpp \users\jang\c\lib\wave\waveRead4pda.cpp -I/users/jang/c/lib -I/users/jang/c/lib/audio -I/users/jang/c/lib/utility -I/users/jang/c/lib/wave

pv=[NaN;48.6873758270755;48.9024387230427;48.6873758270755;48.5808380109181;48.474951806274;48.474951806274;48.5808380109181;48.5808380109181;48.5808380109181;49.0109803982788;49.0109803982788;49.1202068835658;49.0109803982788;48.9024387230427;49.0109803982788;NaN;NaN;NaN;NaN;NaN;NaN;NaN;53.7392170143959;53.8828887143888;54.0277627019957;54.9234283127088;55.3892592993951;55.7069258357294;55.8679720248773;56.1946296496976;56.6964717580331;57.2132964207357;57.3890577323085;57.5666217554455;57.5666217554455;57.7460258530372;58.1105096415594;57.9273085617197;58.1105096415594;58.1105096415594;58.1105096415594;58.1105096415594;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
%pv(isnan(pv))=0;
fs=8000;
pitchRate=100.23;

%pitchPlay(pv, pitchRate, fs, 1);

subplot(2,1,1)
plot(pv, '.-'); grid on; axis tight
wave=pitch2waveMex(pv, pitchRate, fs);
subplot(2,1,2);
plot(wave); axis tight

%sound(wave, fs);

return

% Read a midi file, convert it to PV, synthesize the corresponding wave, and play the wave
% This breaks because readmidi (in midi toolbox) does not support 64-bit WIN7 platform!!!
midiFile='london_bridge.mid';
pitchRate=8000/256;
plotOpt=1;
pv=midiFile2pv(midiFile, pitchRate, plotOpt);
wave=pitch2waveMex(pv, pitchRate, fs);
sound(wave, fs);
