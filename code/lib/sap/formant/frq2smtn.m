function smtn = frq2smtn(freq)
% FRQ2SMTN Frequency to MIDI semitone conversion
smtn = 69+12*log2(freq/440);
