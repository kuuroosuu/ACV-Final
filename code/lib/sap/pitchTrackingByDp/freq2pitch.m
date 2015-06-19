function pitch=freq2pitch(freq)
% freq2pitch: Frequency (in Hz) to pitch (in semitone) conversion

%	Roger Jang, 20040524

pitch = 69+12*log2(freq/440);