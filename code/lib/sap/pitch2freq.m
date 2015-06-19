function freq = pitch2freq(pitch)
% pitch2freq: Pitch (in semitones) to freqency (in Hertz) conversion. [由音高（以Semitone為單位）至頻率（以Hertz為單位）的轉換]
%	Usage: freq = pitch2freq(pitch)
%		pitch: pitch vector
%		freq: frequency vector
%		(Elements of zero are considered as silence and thus not converted.)

%	Roger Jang, 2004xxxx

freq=pitch;
index=find(pitch~=0);
freq(index) = 440*2.^((pitch(index)-69)/12);
