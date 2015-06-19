function notePlot(note)
%notePlot: Monophonic music note plotting
%	Usage: notePlot(note)
%		note: music note vector with two fields
%			note.pitch: pitch in semitones, where zero/nan is considered as silence
%			note.duration: duration in seconds

%	Roger Jang, 20060303

if nargin<1, selfdemo; return; end

pitch=[note.pitch];
pitch=pitch(:)';
duration=[note.duration];
duration=duration(:)';

%pitch(pitch==0)=nan;
plot([0, cumsum(duration)], [pitch, nan], 'r.');
cumDuration = [0 cumsum(duration)];
tmp = [cumDuration(1:end-1); cumDuration(2:end)];
t = tmp(:);
tmp = [pitch; pitch];
x = tmp(:);
lineH = line(t, x, 'color', 'b', 'linestyle', '-');
xlabel('Duration (Seconds)');
ylabel('Pitch (Semitones)');
grid on
axis([-inf, inf, min(pitch)-1, max(pitch)+1]);

% ====== Self demo
function selfdemo
pitch=   [55 55 55 55 57 55 0 57 60 0  64 0  62 62 62 62 60 64 60 0 57 55 55];
duration=[23 23 23 23 23 35 9 23 69 18 69 18 23 23 23 12 12 23 35 9 12 12 127]/64;
for i=1:length(pitch)
	note(i).pitch=pitch(i);
	note(i).duration=duration(i);
end
feval(mfilename, note);
title('「哭砂」前兩句的音符');
