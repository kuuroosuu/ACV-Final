note(1).pitch=60; note(1).duration=0.3;
note(2).pitch=62; note(2).duration=0.3;
note(3).pitch=64; note(3).duration=0.3;

pitch=   [55 55 55 55 57 55 0 57 60 0  64 0  62 62 62 62 60 64 60 0 57 55 55];
duration=[23 23 23 23 23 35 9 23 69 18 69 18 23 23 23 12 12 23 35 9 12 12 127]/64;
for i=1:length(pitch)
	note(i).pitch=pitch(i);
	note(i).duration=duration(i);
end

pitch=[note.pitch];
duration=[note.duration];
fs=16000;
wave=note2waveMex(pitch, duration, 16000);
%sound(wave, fs);
close all
plot(wave);