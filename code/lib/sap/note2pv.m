function [pv, noteStartPvIndex, anchorPvIndex]=note2pv(note, pitchRate, pvTotalDuration, anchorNoteIndex, plotOpt)
%note2pv: note to pv conversion
%	Usage: [pv, noteStartPvIndex, anchorPvIndex] = note2pv(note, pitchRate, pvTotalDuration, anchorNoteIndex, plotOpt)
%		note: note vector of [pitch1, duration1, pitch2, duration2, ...], where the unit of duration is 1/64 seconds
%		pitchRate: pitch points per second
%		pvTotalDuration: total duration of the returned PV (inf for the whole note vector)
%		anchorNoteIndex: indices of notes of anchored position (where the match begins)
%		pv: pitch vector
%
%	Type note2pv for a self demo.

%	Roger Jang, 20040612, 20090914

if nargin<1, selfdemo; return; end
if nargin<2, pitchRate=8000/256; end
if nargin<3, pvTotalDuration=inf; end
if nargin<4, anchorNoteIndex=[]; end
if nargin<5, plotOpt=0; end

note=double(note);
pitch=note(1:2:end);		% Semitone
duration=note(2:2:end)/64;	% Duration in second

noteIndex=1;
cumDuration=duration(noteIndex);
pv=pitch(noteIndex);
t=0;
count=1;
noteStartPvIndex=[1];		% Index (in pv) of note start positions
anchorPvIndex=[1];		% Index (in pv) of anchor positions
pvTimeStep=1/pitchRate;
while 1
	t=t+pvTimeStep;
	count=count+1;
	if t>pvTotalDuration
		break;
	end
	if t>=cumDuration
		noteIndex=noteIndex+1;
		if noteIndex>length(pitch)
			break;
		end
		noteStartPvIndex=[noteStartPvIndex, count];
		if ~isempty(find(anchorNoteIndex==noteIndex))
			anchorPvIndex=[anchorPvIndex, count];
		end
		cumDuration=cumDuration+duration(noteIndex);
	end
	pv=[pv, pitch(noteIndex)];
end
% Only keep noteStartPvIndex which generate non-zero semitone value
pvAtNoteStart=pv(noteStartPvIndex);
noteStartPvIndex(pvAtNoteStart==0)=[];

if plotOpt
	temp=pv;
	temp(temp==0)=nan;
	pvTime=pvTimeStep*(0:length(temp)-1);
	plot(pvTime, temp, 'g.');
	hold on
	notePlot(note);
	hold off
	% Use a triangle to indicate the starting pv of a note
	for i=1:length(noteStartPvIndex)
		line(pvTime(noteStartPvIndex(i)), pv(noteStartPvIndex(i)), 'marker', '^', 'color', 'k');
	end
	% Use a bigger triangle to indicate the starting pv of a anchor note
	for i=1:length(anchorPvIndex)
		line(pvTime(anchorPvIndex(i)), pv(anchorPvIndex(i)), 'marker', '^', 'color', 'r', 'markersize', 12);
	end
end

% ====== Self demo
function selfdemo
% 南海姑娘
pitch=[55 60 64 62 60 55 57 0 55 60 65 64 62 60 67 67 69 72 71 67 71 69 65 69 67 64 65 67 65 60 64 62 55 62 60];
duration=0.5*[0.75 0.25 0.5 0.5 0.5 0.5 4 1 0.5 0.5 0.5 0.5 0.5 0.5 3 0.5 0.5 0.5 1 0.5 0.5 0.5 0.5 0.5 3 0.5 0.5 0.5 1 0.5 0.5 0.5 0.5 0.5 4]*64;
note=[pitch; duration]; note=note(:)';
% 「哭砂」的前兩句
%note = [55 23 55 23 55 23 55 23 57 23 55 35 0 9 57 23 60 69 0 18 64 69 0 18 62 23 62 23 62 23 62 12 60 12 64 23 60 35 0 9 57 12 55 12 55 127];
pitchRate=256/8000;
pvTotalDuration=9;
anchorNoteIndex=[1 5 12];
[pv, noteStartPvIndex, anchorPvIndex]=note2pv(note, pitchRate, pvTotalDuration, anchorNoteIndex, 1);
fprintf('noteStartPvIndex=%s\n', mat2str(noteStartPvIndex));
fprintf('anchorPvIndex=%s\n', mat2str(anchorPvIndex));
title(['Green dots is the output ', mfilename]);
fprintf('Press any key to play note...\n'); pause
notePlay(note);
fprintf('Press any key to play pv...\n'); pause
pvPlay(pv, pitchRate);