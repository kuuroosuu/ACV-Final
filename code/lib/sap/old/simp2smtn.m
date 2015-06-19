function smtn = simp2smtn(note, octave)
% NOTE2SEM Conversion from simplified note representation (Â²ÃÐ) to semitone.

%	Roger Jang, 990928

if nargin==0, selfdemo; return; end

note_n = length(note);
simple2semitone = [-9 -7 -5 -4 -2 0 2];
smtn = zeros(1, note_n);
for i = 1:note_n,
	if note(i)==0,
		smtn(i) = nan;
	elseif note(i)~=floor(note(i)),
		below = floor(note(i));
		above = ceil(note(i));
		smtn(i) = mean(simple2semitone([above, below]));
	elseif (1<=note(i)) & (note(i)<=7),
		smtn(i) = simple2semitone(note(i));
	else
		error('Error in processing the note!');
	end
end

% Modify semitones according to octave
smtn = smtn + 12*octave + 69;

% ========= Selfdemo ==========
function selfdemo
simple = [0 5 1 3 2 1 5 6 0 5 1 4 3 2 1 5 5 6 1 7 5 7 6 4 6 5 3 4 5 4 1 3 2 5 2 1];
octave = [0 -1 0 0 0 0 -1 -1 0 -1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0];
beat = [1 .75 .25 .5 .5 .5 .5 4 1 .5 .5 .5 .5 .5 .5 3 .5 .5 .5 1 .5 .5 .5 .5 .5 3 .5 .5 .5 1 .5 .5 .5 .5 .5 4]*500;

smtn = feval(mfilename, simple, octave);
note = reshape([smtn; beat], 1, 2*length(beat));
fs = 11025;
% ====== Plot the wave file
wave = midi2wave(note, fs);
subplot(2,1,1); plot((1:length(wave))/fs, wave);
xlabel('Time (sec)');
ylabel('Wave intensity');
axis tight;
% ====== Plot the pitch contour
subplot(2,1,2);
plotmidi(note, 1/1000);
% ====== Play the sound wave
fprintf('Hit return to play...\n');
pause
sound(wave, fs);
%playmidi(note, 1/1000);
