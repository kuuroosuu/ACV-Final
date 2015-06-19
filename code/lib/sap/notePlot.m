function lineH = notePlot(note, timeUnit, color, startTime)
%notePlot: Plot note-format music sequence 
%	Usage: notePlot(note, timeUnit)
%		note(1:2:end): Semitone in MIDI representation
%		note(2:2:end): Duration in timeUnit second
%		(Default value for timeUnit is 1/64 seconds.)
%
%	Example:
%		pitch=[68 65 65 66 63 63 61 63 65 66 68 68 68];
%		duration=[1 1 2 1 1 2 1 1 1 1 1 1 2]*0.5*64;
%		note=reshape([pitch; duration], 1, 2*length(pitch));
%		notePlot(note, 1/64);
%		title('Notes for "Little Bee"');
%		notePlay(note, 1/64);
%
%	You can try "notePlot" to see the self demo.

%	Category: Visualization
%	Roger Jang, 19990928, 20010204, 20070601

if nargin<1, selfdemo; return; end
if nargin<2, timeUnit=1/64; end
if nargin<3, color='b'; end
if nargin<4, startTime=0; end

if isstruct(note)
	pitch=[note.pitch];
	duration=[note.duration]*timeUnit;
else
	pitch=double(note(1:2:end));
	duration=double(note(2:2:end))*timeUnit;
end
pitch(pitch==0)=nan;

pitch=pitch(:)';
duration=duration(:)';

cumtime = [0 cumsum(duration)];
tmp = [cumtime(1:end-1); cumtime(2:end)];
t = tmp(:)+startTime;
tmp = [pitch; pitch];
x = tmp(:);
lineH = plot(t, x);
line([0, cumsum(duration)]+startTime, [pitch, nan], 'color', 'r', 'marker', 'o', 'linestyle', 'none');
xlabel('Time (Seconds)');
ylabel('Pitch (Semitones)');
grid on
%axis tight;

% ========= Selfdemo ==========
function selfdemo
pitch=[0 55 60 64 62 60 55 57 0 55 60 65 64 62 60 67 67 69 72 71 67 71 69 65 69 67 64 65 67 65 60 64 62 55 62 60];
duration=[1 0.75 0.25 0.5 0.5 0.5 0.5 4 1 0.5 0.5 0.5 0.5 0.5 0.5 3 0.5 0.5 0.5 1 0.5 0.5 0.5 0.5 0.5 3 0.5 0.5 0.5 1 0.5 0.5 0.5 0.5 0.5 4]*0.5*64;
note=[pitch; duration]; note=note(:)';
notePlot(note, 1/64);
title('Music notes for "Girl Over Southern Ocean" (南海姑娘)');
fprintf('Hit return to play "Girl Over Southern Ocean" (南海姑娘):\n'); pause
notePlay(note, 1/64);